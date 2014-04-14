//(C) Microsoft Corporation. All rights reserved.
#include "stdafx.h"
#include "Graph.h"
#include "Utils.h"
#include "mex.h"

namespace paraF
{

const double GLOBAL_RELABEL_ARC_COEFF = 1;
const double GLOBAL_RELABEL_NODE_COEFF = 6;

const int RELABEL_WORK_CONST = 12;
const int RELABEL_WORK_PER_ARC = 1;

const int GLOBAL_UPDATE_WORK_CONST   = 5;
const int GLOBAL_UPDATE_WORK_PER_ARC = 2;

const int PUSH_WORK_CONST = 0;

// Initializes the graph. Sets multiplier to 1. Resets counters.
Graph::Graph()
{
	multiplier = 1;
	ClearStats();
}

// Clears statistics.
void Graph::ClearStats()
{
	pushCounter = 0;
	relabelCounter = 0;
	gapCounter = 0;
	globalUpdateCounter = 0;
	cleanupFetchCounter = 0;
	workCounter = 0;
	workSinceLastUpdate = 0;
	relabelFetchCounter = 0;
}

// Reserves memory for nodes, arcs etc
void Graph::Reserve(int m, int n)
{
	nodes.reserve(n);
	arcs.reserve(m);
	activeBuckets.resize(n);
	inactiveBuckets.resize(n);
}

// Returns the number of node in the graph (after contractions).
int Graph::GetNodeCount()
{
	return state.firstSinkIndex - state.firstInnerIndex + 2;
}

// Adds a node to the graph.
Node* Graph::AddNode()
{
	if (nodes.size() == nodes.capacity())
		throw "too many nodes";

	Node node;
	node.excess = 0;
	node.height = -1;
	node.capDelta = 0;
	node.svArc = node.vtArc = NULL;
	node.svCap = node.vtCap = node.svSlope = node.vtSlope = 0;
	node.typedIndex = -1;
	node.index = (int) nodes.size();
	nodes.push_back(node);
	return &nodes[node.index];
}

// Adds an arc to the graph.
Arc* Graph::AddArc(Node* tailNode, Node* headNode, int cap)
{
	if (arcs.size() == arcs.capacity())
		throw "too many arcs";

	Arc arc;
	arc.cap = cap;
	arc.headNode = headNode;
#ifdef VERBOSE
	arc.index = (int) arcs.size();
#endif
	arcs.push_back(arc);
	return &arcs[arcs.size() - 1];
}

// Adds an arc leaving source node to the graph.
Arc* Graph::AddArcFromSource(Node* tailNode, Node* headNode, int cap, int slope)
{
	Arc* arc = AddArc(tailNode, headNode, cap);
	if (headNode->svArc != NULL)
		throw "parallel arcs from source are disallowed";
	headNode->svArc = arc;
	headNode->svCap = cap;
	headNode->svSlope = slope;
	return arc;
}

// Adds an arc entring sink node to the graph.
Arc* Graph::AddArcToSink(Node* tailNode, Node* headNode, int cap, int slope)
{
	Arc* arc = AddArc(tailNode, headNode, cap);
	if (tailNode->vtArc != NULL)
		throw "parallel arcs to sink are disallowed";
	tailNode->vtArc = arc;
	tailNode->vtCap = cap;
	tailNode->vtSlope = slope;
	return arc;
}

// Compares two arcs by tail node index.
class ArcComparer
{
public:
	ArcComparer(Graph* graph)
	{
		this->graph = graph;
	}

	bool operator () (int lhs, int rhs) const
	{
		Arc& lhsArc = graph->arcs[lhs];
		Arc& rhsArc = graph->arcs[rhs];
		return lhsArc.GetTailNode()->index < rhsArc.GetTailNode()->index;
	}

private:
	Graph* graph;

};

// Sorts arcs by tail node, prepares outArcBegin/outArcEnd.
void Graph::PrepareArcs()
{
	int m = (int) arcs.size();
	std::vector<int> perm(m);
	std::vector<int> invPerm(m);
	for (int i = 0; i < m; i++)
		perm[i] = i;

	std::sort(perm.begin(), perm.end(), ArcComparer(this));

	for (int i = 0; i < m; i++)
		invPerm[perm[i]] = i;

	for (int i = 0; i < m; i++)
		if (perm[i] >= 0)
		{
			Arc curArc = arcs[i];
			arcs[i] = arcs[perm[i]];
			perm[i] = -1;
			int j = invPerm[i];
			while (perm[j] >= 0)
			{
				Arc newCurArc = arcs[j];
				arcs[j] = curArc;
				perm[j] = -1;
				curArc = newCurArc;
				j = invPerm[j];
			}
		}

	for (ArcIterator arc = arcs.begin(); arc != arcs.end(); arc++)
	{
		arc->revArc = &arcs[invPerm[arc->revArc - &*arcs.begin()]];
	}

	for (NodeIterator node = nodes.begin(); node != nodes.end(); node++)
	{
		node->svArc = node->svArc == NULL ? NULL : &arcs[invPerm[node->svArc - &*arcs.begin()]];
		node->vtArc = node->vtArc == NULL ? NULL : &arcs[invPerm[node->vtArc - &*arcs.begin()]];
		node->outArcBegin = NULL;
		node->outArcEnd = NULL;
	}

	Node* curNode = &*nodes.begin();
	for (ArcIterator arc = arcs.begin(); arc != arcs.end(); arc++)
	{
		if (arc->GetTailNode() != curNode)
		{
			curNode->outArcEnd = &*arc;
			curNode++;
		}

		if (arc->GetTailNode() == curNode && curNode->outArcBegin == NULL)
			curNode->outArcBegin = &*arc;
	}

	if (curNode->outArcBegin != NULL)
		curNode->outArcEnd = &*arcs.begin() + arcs.size();
}

// Initializes nodes and typedNodeList.
void Graph::PrepareNodes(Node* source, Node* sink)
{
	for (NodeIterator node = nodes.begin(); node != nodes.end(); node++)
	{
		node->typedIndex = (int) typedNodeList.size();
		typedNodeList.push_back(&*node);
	}

	state.firstInnerIndex = 1;
	Node* node0 = typedNodeList[0];
	std::swap(typedNodeList[0], typedNodeList[source->typedIndex]);
	node0->typedIndex = source->typedIndex;
	source->typedIndex = 0;

	state.firstSinkIndex = (int) nodes.size() - 1;
	Node* node1 = typedNodeList[state.firstSinkIndex];
	std::swap(typedNodeList[state.firstSinkIndex], typedNodeList[sink->typedIndex]);
	node1->typedIndex = sink->typedIndex;
	sink->typedIndex = state.firstSinkIndex;

	state.sourceCap = 0;
	state.sourceSlope = 0;
	state.sinkCap = 0;
	state.sinkSlope = 0;
	for (int i = state.firstInnerIndex; i < state.firstSinkIndex; i++)
	{
		Node* node = typedNodeList[i];
		state.sourceCap += node->svCap;
		state.sourceSlope += node->svSlope;
		state.sinkCap += node->vtCap;
		state.sinkSlope += node->vtSlope;
	}
}

// Checks if a given node is one of the source nodes (in the contracted graph).
bool Graph::IsSourceNode(Node* node)
{
	return node->typedIndex < state.firstInnerIndex;
}

// Checks if a given node is one of the sink nodes (in the contracted graph).
bool Graph::IsInnerNode(Node* node)
{
	return node->typedIndex >= state.firstInnerIndex && node->typedIndex < state.firstSinkIndex;
}

// Checks if a given node is one of the inner nodes (in the contracted graph).
bool Graph::IsSinkNode(Node* node)
{
	return node->typedIndex >= state.firstSinkIndex;
}

// Constructs maximum preflow in the graph or exists when
// more than maxWorkAmount units of time is spent.
bool Graph::FindMaxPreflow(int maxWorkAmount)
{
	while (workCounter < maxWorkAmount)
	{
		if (workSinceLastUpdate > globalUpdateBarrier)
		{
			GlobalUpdate();
			workSinceLastUpdate -= globalUpdateBarrier;
		}

		Node* node = PopNodeForDischarge();
		if (node == NULL)
			return true;

		LongType workCounterBefore = workCounter;
		DischargeNode(node);
		workSinceLastUpdate += (workCounter - workCounterBefore);
	}
	workCounter -= maxWorkAmount;
	return false;
}

// Constructs maximum preflow (never exists prematurely).
void Graph::FindMaxPreflow()
{
	FindMaxPreflow(std::numeric_limits<int>::max());
}

// Pushes a given amount of preflow through a given arc.
void Graph::PushPreflow(Arc* arc, LongType amount)
{
#ifdef VERBOSE
	printf("c Pushing %" LONGTYPE_SPEC " units through arc %d->%d\n", amount, arc->GetTailNode()->index, arc->GetHeadNode()->index);
#endif

	workCounter += PUSH_WORK_CONST;
	pushCounter++;

	arc->flow += amount;
	arc->revArc->flow -= amount;

	Node* tailNode = arc->GetTailNode();
	tailNode->excess -= amount;
#ifdef PARANOID
	if (tailNode != GetSourceNode() && tailNode->excess < 0)
		throw "excess underflow";
#endif

	Node* headNode = arc->GetHeadNode();
	if (headNode->excess == 0 && IsInnerNode(headNode))
	{
		RemoveFromInactiveBucket(headNode);
		AddToActiveBucket(headNode);
	}
	headNode->excess += amount;
#ifdef PARANOID
	if (headNode->excess < 0)
		throw "excess overflow";
#endif
}

// Initializes heights. Sets height to n at sources, 1 at inner nodes, and 0 at sinks. 
void Graph::InitHeights()
{
	ClearBuckets();

	GetSourceNode()->height = (int) nodes.size();
	GetSinkNode()->height = 0;
	for (int i = state.firstInnerIndex; i < state.firstSinkIndex; i++)
	{
		Node* node = typedNodeList[i];
		node->height = 1;
		AddToInactiveBucket(&*node);
		node->curArc = node->firstOutArc;
	}
}

// Calculates delta capacities for the graph given the assumption that
// lambda will be greater or equal to a given value.
// JL: when lambda is -1 it seems to just set capDelta to svSlope
void Graph::InitCapDelta(LongType lambda)
{
	for (NodeIterator node = nodes.begin(); node != nodes.end(); node++)
	{
		LongType svCap = node->svCap + lambda * node->svSlope;
		if (svCap < 0)
			node->capDelta = std::max(node->capDelta, -svCap);

		LongType vtCap = node->vtCap + lambda * node->vtSlope;
		if (vtCap < 0)
			node->capDelta = std::max(node->capDelta, -vtCap);
	}
}

// Computes the original capacity of an arc.
// For arcs having parametrization this capacity corresponds to lambda = 0.
LongType Graph::GetArcBaseCap(Arc* arc)
{
	if (arc == arc->GetHeadNode()->svArc)
		return arc->GetHeadNode()->svCap;

	if (arc == arc->GetTailNode()->vtArc)
		return arc->GetTailNode()->vtCap;

	return arc->cap;
}

// Initializes arc capacities for a given lambda taking delta capacities into account.
void Graph::InitCaps(LongType lambda)
{
	for (int i = state.firstInnerIndex; i < state.firstSinkIndex; i++)
	{
		Node* node = typedNodeList[i];
		if (node->svArc != NULL)
		{
			node->svArc->cap = (LongType) node->svCap + lambda * node->svSlope + node->capDelta;
#ifdef PARANOID
			if (node->svArc->cap < 0)
				throw "negative capacity";
#endif
		}
		if (node->vtArc != NULL)
		{
			node->vtArc->cap = (LongType) node->vtCap + lambda * node->vtSlope + node->capDelta;
#ifdef PARANOID
			if (node->vtArc->cap < 0)
				throw "negative capacity";
#endif
		}
	}
}

// Adds a given arc to the list of alive arcs for this tail node.
void Graph::AddAliveArc(Arc* arc)
{
	arc->nextOutArc = arc->GetTailNode()->firstOutArc;
	arc->GetTailNode()->firstOutArc = arc;
}

// Initializes alive arcs list.
void Graph::InitArcLists()
{
	Node* source = GetSourceNode();
	Node* sink = GetSinkNode();

	source->firstOutArc = NULL;
	sink->firstOutArc = NULL;
	for (int i = state.firstInnerIndex; i < state.firstSinkIndex; i++)
	{
		typedNodeList[i]->firstOutArc = NULL;
	}

	for (int i = state.firstInnerIndex; i < state.firstSinkIndex; i++)
	{
		Node* node = typedNodeList[i];
		AddAliveArc(node->svArc);
		AddAliveArc(node->vtArc->revArc);
		for (Arc* arc = node->outArcBegin; arc != node->outArcEnd; arc++)
		{
			cleanupFetchCounter++;
			Node* headNode = arc->GetHeadNode();
			if (IsInnerNode(headNode) || arc == node->vtArc || arc == node->svArc->revArc)
				AddAliveArc(arc);
		}
	}
}

// Clears lists of alive arcs.
void Graph::ClearAliveArcs(Node* node)
{
	for (Arc* arc = node->firstOutArc; arc != NULL; arc = arc->nextOutArc)
	{
		arc->flow = arc->revArc->flow = 0;
	}
}

// Initializes the whole preflow push algorithm using a given value of lambda.
void Graph::InitPreflowPush(LongType lambda)
{
	this->state.lambda = lambda;

	InitArcLists();

	Node* source = GetSourceNode();
	Node* sink = GetSinkNode();
	
	ClearAliveArcs(source);
	ClearAliveArcs(sink);
	for (int i = state.firstInnerIndex; i < state.firstSinkIndex; i++)
	{
		Node* node = typedNodeList[i];
		node->excess = 0;
		ClearAliveArcs(node);
	}

	InitCaps(lambda);
	ClearBuckets();
	InitHeights();

	for (Arc* arc = source->firstOutArc; arc != NULL; arc = arc->nextOutArc)
	{
		LongType resCap = arc->cap - arc->flow;
		if (resCap > 0)
			PushPreflow(arc, resCap);
	}

	globalUpdateBarrier = (int) (GLOBAL_RELABEL_ARC_COEFF * arcs.size() + GLOBAL_RELABEL_NODE_COEFF * nodes.size());

#ifdef VERBOSE
	printf("c Preflow initialized for lambda = %.3lf\n", (double) lambda/multiplier);
#endif
}

// Increases the value of lambda. Adjusts the flow in order to satisfy
// capacity constraints.
void Graph::IncreaseLambda(LongType lambda)
{
#ifdef VERBOSE
	printf("c Increasing lambda to %.3lf\n", (double) lambda/multiplier);
#endif

	state.lambda = lambda;

	InitCaps(lambda);
	InitHeights();

	Node* source = GetSourceNode();
	for (Arc* arc = source->firstOutArc; arc != NULL; arc = arc->nextOutArc)
	{
		arc->flow = arc->cap;
		arc->revArc->flow = -arc->flow;
	}

	for (int i = state.firstInnerIndex; i < state.firstSinkIndex; i++)
	{
		Node* node = typedNodeList[i];
		LongType excess = 0;
		for (Arc* arc = node->firstOutArc; arc != NULL; arc = arc->nextOutArc)
		{
			if (arc != node->vtArc)
				excess -= arc->flow;
		}
#ifdef PARANOID
		if (excess < 0)
			throw "negative excess";
#endif
		LongType amount = std::min(excess, node->vtArc->cap);
		node->vtArc->flow = amount;
		node->vtArc->revArc->flow = -amount;
		excess -= amount;
		node->excess = excess;
	}

	Node* sink = GetSinkNode();
	for (Arc* arc = sink->firstOutArc; arc != NULL; arc = arc->nextOutArc)
	{
		if (IsInnerNode(arc->GetHeadNode()))
		{
			Arc* revArc = arc->revArc;
			LongType amount = revArc->flow - revArc->cap;
			if (amount > 0)
				PushPreflow(arc, amount);
		}
	}

	GlobalUpdate();
}

// Computes a lower estimate for lambda.
LongType Graph::GetLowerLambdaBound()
{
	LongType minNom = std::numeric_limits<int>::max();
	LongType minDenom = 1;

	for (int i = state.firstInnerIndex; i < state.firstSinkIndex; i++)
	{
		Node* node = typedNodeList[i];

		LongType denom = node->svSlope - node->vtSlope;
		if (denom <= 0)
			continue;
		
		LongType nom = node->vtCap - node->svCap;
		for (Arc* arc = node->outArcBegin; arc != node->outArcEnd; arc++)
		{
			if (!IsInnerNode(arc->GetHeadNode()))
				continue;
			Arc* revArc = arc->revArc;
			nom -= revArc->cap / multiplier;
		}

		LongType g = gcd(abs(nom), abs(denom));
		nom /= g;
		denom /= g;

		if (nom * minDenom - denom * minNom < 0)
		{
			minNom = nom;
			minDenom = denom;
		}
	}
	return multiplier * minNom / minDenom - multiplier;
}

// Computes an upper estimate for lambda.
LongType Graph::GetUpperLambdaBound()
{
	LongType maxNom = std::numeric_limits<int>::min();
	LongType maxDenom = 1;

	for (int i = state.firstInnerIndex; i < state.firstSinkIndex; i++)
	{
		Node* node = typedNodeList[i];

		LongType denom = node->svSlope - node->vtSlope;
		if (denom <= 0)
			continue;
		
		LongType nom = node->vtCap - node->svCap;
		for (Arc* arc = node->outArcBegin; arc != node->outArcEnd; arc++)
		{
			if (!IsInnerNode(arc->GetHeadNode()))
				continue;
			nom += arc->cap / multiplier;
		}

		LongType g = gcd(abs(nom), abs(denom));
		nom /= g;
		denom /= g;

		if (nom * maxDenom - denom * maxNom > 0)
		{
			maxNom = nom;
			maxDenom = denom;
		}
	}
	return multiplier * maxNom / maxDenom + multiplier;
}

// Discharges a given node.
void Graph::DischargeNode(Node* node)
{
#ifdef VERBOSE
 	printf("c Discharing node %d\n", node->index);
#endif

	for (;;)
	{
		Arc* curArc = node->curArc;
		LongType resCap = curArc->cap - curArc->flow;
		if (resCap > 0)
		{
			Node* headNode = curArc->GetHeadNode();
			if (headNode->height < node->height)
			{
				LongType amount = std::min(node->excess, resCap);
				PushPreflow(curArc, amount);
				if (node->excess == 0)
					break;
			}
		}
		node->curArc = curArc->nextOutArc;
		if (node->curArc == NULL)
		{
			if (!RelabelNode(node))
				break;
		}
	}

	if (node->height < (int) nodes.size())
	{
		if (node->excess > 0)
			AddToActiveBucket(node);
		else
			AddToInactiveBucket(node);
	}
}

// Relabels a given node.
bool Graph::RelabelNode(Node* node)
{
	relabelCounter++;
	workCounter += RELABEL_WORK_CONST;

	int oldHeight = node->height;
	if (activeBuckets[oldHeight].firstNode == NULL && inactiveBuckets[oldHeight].firstNode == NULL)
	{
		Gap(node->height);
		node->height = (int) nodes.size();
		return false;
	}

	int minHeight = std::numeric_limits<int>::max();
	Arc* minArc;
	for (Arc* arc = node->firstOutArc; arc != NULL; arc = arc->nextOutArc)
	{
		workCounter += RELABEL_WORK_PER_ARC;
		relabelFetchCounter++;
		int curHeight = arc->GetHeadNode()->height;
		LongType resCap = arc->cap - arc->flow;
		if (resCap > 0 && minHeight > curHeight)
		{
			minHeight = curHeight;
			minArc = arc;
		}
	}

	minHeight++;
	node->curArc = minArc;

#ifdef VERBOSE
	printf("c Relabeled node %d from %d to %d\n", node->index, node->height, minHeight);
#endif

	node->height = minHeight;
	return minHeight < (int) nodes.size();
}

// Adds a given node to the list of active nodes (corresponding to its bucket).
void Graph::AddToActiveBucket(Node* node)
{
	int height = node->height;
	Bucket& bucket = activeBuckets[height];
	node->nextInBucket = bucket.firstNode;
	bucket.firstNode = &*node;
	maxHeight = std::max(maxHeight, height);
}

// Adds a given node to the list of inactive nodes (corresponding to its bucket).
void Graph::AddToInactiveBucket(Node* node)
{
	int height = node->height;
	Bucket& bucket = inactiveBuckets[height];
	if (bucket.firstNode == NULL)
	{
		node->prevInBucket = node->nextInBucket = NULL;
		bucket.firstNode = &*node;
	}
	else
	{
		node->nextInBucket = bucket.firstNode;
		node->prevInBucket = NULL;
		bucket.firstNode->prevInBucket = &*node;
		bucket.firstNode = &*node;
	}
}

// Removes a given node from the list of inactive nodes (corresponding to its bucket).
void Graph::RemoveFromInactiveBucket(Node* node)
{
	int height = node->height;
	Bucket& bucket = inactiveBuckets[height];
	if (node->prevInBucket == NULL)
		bucket.firstNode = node->nextInBucket;
	else
		node->prevInBucket->nextInBucket = node->nextInBucket;

	if (node->nextInBucket != NULL)
		node->nextInBucket->prevInBucket = node->prevInBucket;
}

// Gets a node for discharge (should be one of active nodes with maximum height).
Node* Graph::PopNodeForDischarge()
{
	while (maxHeight >= 0 && activeBuckets[maxHeight].firstNode == NULL)
		maxHeight--;

	if (maxHeight < 0)
		return NULL;

	Bucket& bucket = activeBuckets[maxHeight];
	Node* node = bucket.firstNode;
	bucket.firstNode = node->nextInBucket;

	return node;
}

// Applies gap relabelling from a given height.
void Graph::Gap(int height)
{
#ifdef VERBOSE
	printf("c Found gap at level %d\n", height);
#endif

	gapCounter++;
	int n = (int) nodes.size();

	for (int i = height; i <= maxHeight; i++)
	{
		Bucket& activeBucket = activeBuckets[i];
		for (Node* node = activeBucket.firstNode; node != NULL; node = node->nextInBucket)
			node->height = n;
		activeBucket.firstNode = NULL;

		Bucket& inactiveBucket = inactiveBuckets[i];
		for (Node* node = inactiveBucket.firstNode; node != NULL; node = node->nextInBucket)
			node->height = n;
		inactiveBucket.firstNode = NULL;
	}
	maxHeight = height - 1;
}

// Clears all buckets.
void Graph::ClearBuckets()
{
	int n = GetNodeCount();
	for (int i = 0; i < n; i++)
	{
		activeBuckets[i].firstNode = NULL;
		inactiveBuckets[i].firstNode = NULL;
	}
	maxHeight = 0;
}

// Updates all heights by performing the backward BFS from the sink node.
void Graph::GlobalUpdate()
{
#ifdef VERBOSE
	printf("c Performing global update\n");
#endif
	globalUpdateCounter++;
	workCounter += GLOBAL_UPDATE_WORK_CONST;

	std::queue<Node*> bfsQueue;
	std::vector<bool> visited(nodes.size());

	Node* sink = GetSinkNode();
	bfsQueue.push(sink);
	visited[sink->index] = true;

	ClearBuckets();

	while (!bfsQueue.empty())
	{
		Node* node = bfsQueue.front();
		bfsQueue.pop();
		int newHeight = node->height + 1;
		for (Arc* arc = node->firstOutArc; arc != NULL; arc = arc->nextOutArc)
		{
			workCounter += GLOBAL_UPDATE_WORK_PER_ARC;
			Arc* revArc = arc->revArc;
			LongType resCap = revArc->cap - revArc->flow;
			if (resCap > 0)
			{
				Node* anotherNode = arc->GetHeadNode();
				if (!visited[anotherNode->index])
				{
					visited[anotherNode->index] = true;
					anotherNode->height = newHeight;
					if (anotherNode->excess > 0)
						AddToActiveBucket(anotherNode);
					else
						AddToInactiveBucket(anotherNode);
					bfsQueue.push(anotherNode);
				}
			}
		}
	}

	int n = (int) nodes.size();
	for (int i = state.firstInnerIndex; i < state.firstSinkIndex; i++)
	{
		Node* node = typedNodeList[i];
		node->curArc = node->firstOutArc;
		if (!visited[node->index])
			node->height = n;
	}
}

// Constructs the reversed graph.
void Graph::Reverse(Graph& revGraph)
{
	if (revGraph.multiplier != 1)
		throw "multiplier is already set";

	int m = (int) arcs.size();
	int n = (int) nodes.size();
	revGraph.Reserve(m, n);

	for (NodeIterator i = nodes.begin(); i != nodes.end(); i++)
		revGraph.AddNode();

	for (ArcIterator arc = arcs.begin(); arc != arcs.end(); arc++)
	{
		revGraph.AddArc(
			&revGraph.nodes[arc->GetTailNode() - &*nodes.begin()],
			&revGraph.nodes[arc->GetHeadNode() - &*nodes.begin()],
			(int) arc->revArc->cap);
	}

	for (int i = 0; i < n; i++)
	{
		revGraph.nodes[i].svCap = nodes[i].vtCap;
		revGraph.nodes[i].svSlope = -nodes[i].vtSlope;
		revGraph.nodes[i].vtCap = nodes[i].svCap;
		revGraph.nodes[i].vtSlope = -nodes[i].svSlope;
	}

	for (int i = 0; i < m; i++)
	{
		Arc& arc = arcs[i];
		Arc& revArc = revGraph.arcs[i];
		revArc.revArc = &revGraph.arcs[arc.revArc - &*arcs.begin()];
	}

	for (int i = 0; i < n; i++)
	{
		revGraph.nodes[i].outArcBegin = &*revGraph.arcs.begin() + (nodes[i].outArcBegin - &*arcs.begin());
		revGraph.nodes[i].outArcEnd = &*revGraph.arcs.begin() + (nodes[i].outArcEnd - &*arcs.begin());

		revGraph.nodes[i].svArc = nodes[i].vtArc == NULL ? NULL : revGraph.arcs[nodes[i].vtArc - &*arcs.begin()].revArc;
		revGraph.nodes[i].vtArc = nodes[i].svArc == NULL ? NULL : revGraph.arcs[nodes[i].svArc - &*arcs.begin()].revArc;
	}

	revGraph.PrepareNodes(&revGraph.nodes[GetSinkNode() - &*nodes.begin()], &revGraph.nodes[GetSourceNode() - &*nodes.begin()]);
}

// Sets multiplier and applies it arc capacities.
void Graph::SetMultiplier(int multiplier)
{
	this->multiplier = multiplier;
	
	for (ArcIterator arc = arcs.begin(); arc != arcs.end(); arc++)
	{
		arc->cap *= multiplier;
	}

	for (NodeIterator node = nodes.begin(); node != nodes.end(); node++)
	{
		node->capDelta *= multiplier;
		node->svCap *= multiplier;
		node->vtCap *= multiplier;
	}

	state.sourceCap *= multiplier;
	state.sinkCap *= multiplier;
}

// Returns the source node (before any contractions).
Node* Graph::GetSourceNode()
{
	return typedNodeList[0];
}

// Returns the sink node (before any contractions).
Node* Graph::GetSinkNode()
{
	return typedNodeList[(int) nodes.size() - 1];
}

// Adds s->v and v->t arcs.
void Graph::AddAuxArcs(Node* source, Node* sink)
{
	for (NodeIterator nodeIt = nodes.begin(); nodeIt != nodes.end(); nodeIt++)
	{
		Node* node = &*nodeIt;

		if (node == source || node == sink)
			continue;

		if (node->svArc == NULL)
		{
			Arc* arc = AddArcFromSource(source, &*node, 0, 0);
			Arc* revArc = AddArc(&*node, source, 0);
			arc->revArc = revArc;
			revArc->revArc = arc;
			node->svArc = arc;
		}
		if (node->vtArc == NULL)
		{
			Arc* arc = AddArcToSink(&*node, sink, 0, 0);
			Arc* revArc = AddArc(sink, &*node, 0);
			arc->revArc = revArc;
			revArc->revArc = arc;
			node->vtArc = arc;
		}
	}
}

// Pushes excesses back to the source.
void Graph::ConvertPreflowToFlow()
{
	GetSourceNode()->color = Node::WHITE_NODE;
	GetSinkNode()->color = Node::WHITE_NODE;
	for (int i = state.firstInnerIndex; i < state.firstSinkIndex; i++)
	{
		Node* node = typedNodeList[i];
		node->color = Node::WHITE_NODE;
	}

	std::vector<Node*> sortedNodes;
	sortedNodes.reserve(nodes.size());

	for (int i = state.firstInnerIndex; i < state.firstSinkIndex; i++)
	{
		Node* node = typedNodeList[i];
		if (node->color == Node::BLACK_NODE || node->excess == 0)
			continue;

		node->parentArc = NULL;
		node->color = Node::GREY_NODE;
		node->curArc = node->firstOutArc;
		DFS(node, sortedNodes);
	}

	for (std::vector<Node*>::reverse_iterator i = sortedNodes.rbegin(); i != sortedNodes.rend(); i++)
	{
		Node* node = *i;
		if (!IsInnerNode(node) || node->excess == 0)
			continue;

		for (Arc* arc = node->firstOutArc; arc != NULL; arc = arc->nextOutArc)
		{
			if (arc->flow > 0)
				continue;
			
			Arc* revArc = arc->revArc;
			LongType amount = std::min(node->excess, -arc->flow);
			arc->flow += amount;
			arc->revArc->flow -= amount;
			node->excess -= amount;
			arc->GetHeadNode()->excess += amount;
		}

		if (node->excess > 0)
			throw "failed to reroute excesses during the second stage";
	}
}

// Auxiliary method used by ConvertPreflowToFlow -- performs depth-first search.
void Graph::DFS(Node* v, std::vector<Node*>& sortedNodes)
{
	for (;;)
	{
		while (!IsSourceNode(v) && v->curArc != NULL)
		{
			Arc* arc = v->curArc;
			Arc* revArc = arc->revArc;
			if (arc->flow >= 0)
			{
				v->curArc = arc->nextOutArc;
				continue;
			}
			
			Node* w = arc->GetHeadNode();
			if (w->color == Node::BLACK_NODE)
			{
				v->curArc = arc->nextOutArc;
				continue;
			}

			if (w->color == Node::GREY_NODE)
			{
				LongType minCap = arc->cap - arc->flow;
				for (Node* u = v; u != w; u = u->parentArc->GetTailNode())
				{
					minCap = std::min(minCap, u->parentArc->cap - u->parentArc->flow);
				}

				arc->flow += minCap;
				arc->revArc->flow -= minCap;
				Node* restartNode = v;
				for (Node* u = v; u != w; u = u->parentArc->GetTailNode())
				{
					Arc* arc = u->parentArc;
					arc->flow += minCap;
					arc->revArc->flow -= minCap;
					if (arc->flow == arc->cap)
						restartNode = arc->GetTailNode();
				}

				for (Node* u = v; u != restartNode; u = u->parentArc->GetTailNode())
				{
					u->color = Node::WHITE_NODE;
				}

				v = restartNode;
			}
			else
			{
				w->curArc = w->firstOutArc;
				w->color = Node::GREY_NODE;
				w->parentArc = arc;
				v = w;
			}
		}

		v->color = Node::BLACK_NODE;
		sortedNodes.push_back(v);

		Arc* arc = v->parentArc;
		if (arc == NULL)
			break;
		v = arc->GetTailNode();
	}
}

// Saves the current state.
void Graph::SaveState(GraphState& state)
{
	state = this->state;
}

// Restores the current state.
void Graph::RestoreState(GraphState& savedState)
{
	while (state.firstInnerIndex > savedState.firstInnerIndex)
	{
		state.firstInnerIndex--;
		Node* node = typedNodeList[state.firstInnerIndex];
		UnmakeSource(node);
	}

	while (state.firstSinkIndex < savedState.firstSinkIndex)
	{
		state.firstSinkIndex++;
		Node* node = typedNodeList[state.firstSinkIndex - 1];
		UnmakeSink(node);
	}

	state = savedState;
}

// Checks if a given arc is alive.
bool Graph::IsArcAlive(Arc* arc)
{
	Node* headNode = arc->GetHeadNode();
	Node* tailNode = arc->GetTailNode();
	Node* sourceNode = GetSourceNode();
	Node* sinkNode = GetSinkNode();

	return
		IsInnerNode(headNode) && IsInnerNode(tailNode) ||
		IsInnerNode(headNode) && arc == headNode->svArc ||
		IsInnerNode(tailNode) && arc == tailNode->vtArc ||
		IsInnerNode(tailNode) && arc->revArc == tailNode->svArc ||
		IsInnerNode(headNode) && arc->revArc == headNode->vtArc;
}

// Adds a given node to the set of sources.
void Graph::MakeSource(Node* node)
{
	state.sourceSlope -= node->svSlope;
	state.sourceSlope += node->vtSlope;
	state.sourceCap -= node->svCap;

	for (Arc* arc = node->outArcBegin; arc != node->outArcEnd; arc++)
	{
		if (!IsArcAlive(arc))
			continue;
		
		Node* headNode = arc->GetHeadNode();
		LongType cap = GetArcBaseCap(arc);
		if (IsInnerNode(headNode))
			headNode->svCap += cap;
		if (!IsSourceNode(headNode))
			state.sourceCap += cap;
	}

	int oldTypedIndex = node->typedIndex;
	if (oldTypedIndex != state.firstInnerIndex)
	{
		std::swap(typedNodeList[state.firstInnerIndex], typedNodeList[oldTypedIndex]);
		node->typedIndex = state.firstInnerIndex;
		typedNodeList[oldTypedIndex]->typedIndex = oldTypedIndex;
	}
	state.firstInnerIndex++;
}

// Restores state changed by preceeding call to MakeSource.
void Graph::UnmakeSource(Node* node)
{
	state.sourceSlope -= node->vtSlope;
	state.sourceSlope += node->svSlope;
	state.sourceCap += node->svCap;

	for (Arc* arc = node->outArcBegin; arc != node->outArcEnd; arc++)
	{
		if (!IsArcAlive(arc))
			continue;
		
		Node* headNode = arc->GetHeadNode();
		LongType cap = GetArcBaseCap(arc);
		if (IsInnerNode(headNode))
			headNode->svCap -= cap;
		if (!IsSourceNode(headNode))
			state.sourceCap -= cap;
	}
}

// Adds a given node to the set of sinks.
void Graph::MakeSink(Node* node)
{
	state.sinkSlope -= node->vtSlope;
	state.sinkSlope += node->svSlope;
	state.sinkCap -= node->vtCap;

	for (Arc* arc = node->outArcBegin; arc != node->outArcEnd; arc++)
	{
		Arc* revArc = arc->revArc;
		if (!IsArcAlive(revArc))
			continue;
		
		Node* tailNode = revArc->GetTailNode();
		LongType cap = GetArcBaseCap(revArc);
		if (IsInnerNode(tailNode))
			tailNode->vtCap += cap;
		if (!IsSinkNode(tailNode))
			state.sinkCap += cap;
	}

	int oldTypedIndex = node->typedIndex;
	if (oldTypedIndex != state.firstSinkIndex - 1)
	{
		std::swap(typedNodeList[state.firstSinkIndex - 1], typedNodeList[oldTypedIndex]);
		node->typedIndex = state.firstSinkIndex - 1;
		typedNodeList[oldTypedIndex]->typedIndex = oldTypedIndex;
	}
	state.firstSinkIndex--;

}

// Restores state changed by preceeding call to MakeSink.
void Graph::UnmakeSink(Node* node)
{
	for (Arc* arc = node->outArcBegin; arc != node->outArcEnd; arc++)
	{
		Arc* revArc = arc->revArc;
		if (!IsArcAlive(revArc))
			continue;
		
		Node* tailNode = revArc->GetTailNode();
		LongType cap = GetArcBaseCap(revArc);
		if (IsInnerNode(tailNode))
			tailNode->vtCap -= cap;
		if (!IsSinkNode(tailNode))
			state.sinkCap -= cap;
	}
}

// Backups current flow values.
void Graph::BackupFlows()
{
	for (int i = state.firstInnerIndex; i < state.firstSinkIndex; i++)
	{
		Node* node = typedNodeList[i];
		for (Arc* arc = node->firstOutArc; arc != NULL; arc = arc->nextOutArc)
		{
			arc->baseFlow = arc->flow;
			arc->revArc->baseFlow = arc->revArc->flow;
		}
	}
}

// Prints performance statistics.
void Graph::PrintStats(char* prefix)
{
	printf("c %s pushes:          %" LONGTYPE_SPEC "\n", prefix, pushCounter);
	printf("c %s relabels:        %" LONGTYPE_SPEC "\n", prefix, relabelCounter);
	printf("c %s gaps:            %" LONGTYPE_SPEC "\n", prefix, gapCounter);
	printf("c %s global updates:  %" LONGTYPE_SPEC "\n", prefix, globalUpdateCounter);
	printf("c %s cleanup fetches: %" LONGTYPE_SPEC "\n", prefix, cleanupFetchCounter);
	printf("c %s relabel fetches: %" LONGTYPE_SPEC "\n", prefix, relabelFetchCounter);
}

// Accumulates performance statistics.
void Graph::AccumulateStats(Graph* graph)
{
	pushCounter += graph->pushCounter;
	relabelCounter += graph->relabelCounter;
	gapCounter += graph->gapCounter;
	globalUpdateCounter += graph->globalUpdateCounter;
	cleanupFetchCounter += graph->cleanupFetchCounter;
	relabelFetchCounter += graph->relabelFetchCounter;
}

// Constructs a miniumum cut. Should be called only when maximum flow is ready.
void Graph::FindMinCut(MinCut& minCut)
{
	minCut.resize(nodes.size());
	for (int i = state.firstInnerIndex; i <  state.firstSinkIndex; i++)
	{
		Node* node = typedNodeList[i];
		minCut[node->index] = Node::INNER_NODE;
	}

	std::queue<Node*> bfsQueue;
	std::vector<bool> visited(nodes.size());

	Node* sink = GetSinkNode();
	bfsQueue.push(sink);
	visited[sink->index] = true;
	minCut[sink->index] = Node::SINK_NODE;

	while (!bfsQueue.empty())
	{
		Node* node = bfsQueue.front();
		bfsQueue.pop();
		for (Arc* arc = node->firstOutArc; arc != NULL; arc = arc->nextOutArc)
		{
			Arc* revArc = arc->revArc;
			LongType resCap = revArc->cap - revArc->flow;
			if (resCap > 0)
			{
				Node* anotherNode = arc->GetHeadNode();
				if (!visited[anotherNode->index])
				{
					visited[anotherNode->index] = true;
					minCut[anotherNode->index] = Node::SINK_NODE;
					bfsQueue.push(anotherNode);
				}
			}
		}
	}

	Node* source = GetSourceNode();
	bfsQueue.push(source);
	visited[source->index] = true;
	minCut[source->index] = Node::SOURCE_NODE;

	while (!bfsQueue.empty())
	{
		Node* node = bfsQueue.front();
		bfsQueue.pop();
		for (Arc* arc = node->firstOutArc; arc != NULL; arc = arc->nextOutArc)
		{
			LongType resCap = arc->cap - arc->flow;
			if (resCap > 0)
			{
				Node* anotherNode = arc->GetHeadNode();
				if (!visited[anotherNode->index])
				{
					visited[anotherNode->index] = true;
					minCut[anotherNode->index] = Node::SOURCE_NODE;
					bfsQueue.push(anotherNode);
				}
			}
		}
	}
}

// Rearranges the nodes in typedNodeList in the following way:
// * the order of sources and sinks is preserved
// * the inner nodes are arranged as follows:
//   - those having cut[i] equal to sourceType,
//   - those having cut[i] equal to Node::INNER_NODE
//   - those having cut[i] equal to -sourceType
// The newFirstInnerIndex parameter is set to the beginning of the second group,
// the newFirstSinkIndex is set to the beginning of the third group.
void Graph::RearrangeNodes(MinCut& cut, Node::MinCutType sourceType, int* newFirstInnerIndex, int* newFirstSinkIndex)
{
	int myFirstInnerIndex = state.firstInnerIndex;
	for (int i = state.firstInnerIndex; i < state.firstSinkIndex; i++)
	{
		Node* node = typedNodeList[i];
		Node::MinCutType type = cut[node - &*nodes.begin()];
		if (type == sourceType)
		{
			int oldTypedIndex = node->typedIndex;
			if (oldTypedIndex != myFirstInnerIndex)
			{
				std::swap(typedNodeList[myFirstInnerIndex], typedNodeList[oldTypedIndex]);
				node->typedIndex = myFirstInnerIndex;
				typedNodeList[oldTypedIndex]->typedIndex = oldTypedIndex;
			}
			myFirstInnerIndex++;
		}
	}
	if (newFirstInnerIndex != NULL)
		*newFirstInnerIndex = myFirstInnerIndex;

	int myFirstSinkIndex = state.firstSinkIndex;
	for (int i = state.firstSinkIndex - 1; i >= myFirstInnerIndex; i--)
	{
		Node* node = typedNodeList[i];
		Node::MinCutType type = cut[node - &*nodes.begin()];
		if (type == (Node::MinCutType) -sourceType)
		{
			myFirstSinkIndex--;
			int oldTypedIndex = node->typedIndex;
			if (oldTypedIndex != myFirstSinkIndex)
			{
				std::swap(typedNodeList[myFirstSinkIndex], typedNodeList[oldTypedIndex]);
				node->typedIndex = myFirstSinkIndex;
				typedNodeList[oldTypedIndex]->typedIndex = oldTypedIndex;
			}
		}
	}
	if (newFirstSinkIndex != NULL)
		*newFirstSinkIndex = myFirstSinkIndex;
}

// Chooses multiplier based on excess estimations.
int Graph::ChooseMultiplier(LongType lambda1, LongType lambda2)
{
	LongType maxEx = 0;
	for (NodeIterator node = nodes.begin(); node != nodes.end(); node++)
	{
		LongType ex = 0;
		for (Arc* arc = node->outArcBegin; arc != node->outArcEnd; arc++)
		{
			Arc* revArc = arc->revArc;
			Node* headNode = revArc->GetHeadNode();
			Node* tailNode = revArc->GetTailNode();
			LongType cap;
			if (revArc == tailNode->vtArc)
			{
				cap =
					std::max((LongType) tailNode->vtSlope * lambda1, (LongType) tailNode->vtSlope * lambda2) +
					tailNode->vtCap + tailNode->capDelta;
			}
			else if (revArc == tailNode->svArc)
			{
				cap =
					std::max((LongType) tailNode->svSlope * lambda1, (LongType) tailNode->svSlope * lambda2) +
					tailNode->svCap + tailNode->capDelta;
			}
			else
				cap = arc->cap;

			if (cap < 0)
				cap = cap;
			ex += cap;
		}
		maxEx = std::max(maxEx, ex);
	}

	return (int) std::min(
		std::numeric_limits<LongType>::max() / 2 / maxEx,
		(LongType) (std::numeric_limits<int>::max() / 2));
}

// Dumps graph information.
void Graph::Dump()
{
	InitArcLists();

	printf("c Lambda = %.3lf, %d node(s), %d arc(s), source nodes = [",\
		(double) state.lambda / multiplier, nodes.size(), arcs.size());
	for (int i = 0; i < state.firstInnerIndex; i++)
		printf(" %d", typedNodeList[i]->index);
	printf(" ], sink nodes = [");
	for (int i = state.firstSinkIndex; i < (int) nodes.size(); i++)
		printf(" %d", typedNodeList[i]->index);
	printf(" ]\n");

	for (NodeIterator nodeIt = nodes.begin(); nodeIt != nodes.end(); nodeIt++)
	{
		Node* node = &*nodeIt;
		if (node != GetSourceNode() && node != GetSinkNode() && !IsInnerNode(node))
			continue;

		printf("c Node %d: height = %d, excess = %" LONGTYPE_SPEC "\n",
			node->index, node->height, node->excess);

		for (Arc* arc = node->firstOutArc; arc != NULL; arc = arc->nextOutArc)
		{
			printf("c  %d -> %d : %.3lf of %.3lf\n",
				arc->GetTailNode()->index, arc->GetHeadNode()->index,
				(double) arc->flow / multiplier, (double) arc->cap / multiplier);
		}
	}

	printf("c Active buckets:\n");
	for (int i = 0; i <= maxHeight; i++)
	{
		Bucket& bucket = activeBuckets[i];
		if (bucket.firstNode != NULL)
		{
			printf("c %d = [", i);
			for (Node* node = bucket.firstNode; node != NULL; node = node->nextInBucket)
				printf(" %d", node->index);
			printf(" ]\n");
		}
	}
	
	printf("c Inactive buckets:\n");
	for (int i = 0; i <= maxHeight; i++)
	{
		Bucket& bucket = inactiveBuckets[i];
		if (bucket.firstNode != NULL)
		{
			printf("c %d = [", i);
			for (Node* node = bucket.firstNode; node != NULL; node = node->nextInBucket)
				printf(" %d", node->index);
			printf(" ]\n");
		}
	}
}

// Checks if preflow is feasible.
void Graph::CheckPreflow()
{
	for (NodeIterator nodeIt = nodes.begin(); nodeIt != nodes.end(); nodeIt++)
	{
		Node* node = &*nodeIt;

		if (node != GetSourceNode() && node != GetSinkNode() && !IsInnerNode(node))
			continue;

		if (IsInnerNode(node) && node->excess < 0)
			throw "negative excess in inner node";

		for (Arc* arc = node->firstOutArc; arc != NULL; arc = arc->nextOutArc)
		{
			if (!IsInnerNode(arc->GetTailNode()) || !IsInnerNode(arc->GetHeadNode()))
				continue;

			if (arc->flow > arc->cap)
				throw "upper flow bounds violated";
			if (arc->flow + arc->revArc->flow != 0)
				throw "antisymmetry violated";
		}
	}
}

// Checks if the flow is feasible.
void Graph::CheckFlow()
{
	CheckPreflow();

	for (int i = state.firstInnerIndex; i < state.firstSinkIndex; i++)
	{
		Node* node = typedNodeList[i];
		LongType sum = 0;
		for (Arc* arc = node->outArcBegin; arc != node->outArcEnd; arc++)
		{
			Node* headNode = arc->GetHeadNode();
			if (IsInnerNode(arc->GetHeadNode()) || headNode == GetSourceNode() || headNode == GetSinkNode())
				sum += arc->flow;
		}
		if (sum != 0)
			throw "conservation violated";
	}
}

// Checks if the flow is feasible and maximum.
void Graph::CheckMaxFlow()
{
	CheckFlow();

	std::queue<Node*> bfsQueue;
	std::vector<bool> visited(nodes.size());

	for (int i = state.firstInnerIndex; i < state.firstSinkIndex; i++)
	{
		Node* node = typedNodeList[i];
		Arc* arc = node->svArc;
		if (arc->flow < arc->cap)
		{
			bfsQueue.push(node);
			visited[node->index] = true;
		}
	}

	while (!bfsQueue.empty())
	{
		Node* node = bfsQueue.front();
		bfsQueue.pop();
		int newHeight = node->height + 1;
		for (Arc* arc = node->firstOutArc; arc != NULL; arc = arc->nextOutArc)
		{
			if (arc->flow < arc->cap)
			{
				Node* anotherNode = arc->GetHeadNode();
				if (!IsSinkNode(anotherNode) && !visited[anotherNode->index])
				{
					visited[anotherNode->index] = true;
					bfsQueue.push(anotherNode);
				}
			}
		}
	}

	for (int i = state.firstInnerIndex; i < state.firstSinkIndex; i++)
	{
		Node* node = typedNodeList[i];
		if (visited[node->index] && node->vtArc->flow < node->vtArc->cap)
			throw "flow is not maximum";
	}
}

}
