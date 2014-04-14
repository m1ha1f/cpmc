//(C) Microsoft corporation. All rights reserved
#ifndef Graph_h_SEEN
#define Graph_h_SEEN


namespace paraF
{

struct Arc;
struct Node;
class Graph;

typedef std::vector<Node>::iterator NodeIterator;
typedef std::vector<Arc>::iterator  ArcIterator;

struct Node
{
	enum NodeColor
	{
		WHITE_NODE,
		GREY_NODE,
		BLACK_NODE
	};

	enum MinCutType
	{
		SOURCE_NODE = 1,
		SINK_NODE = -1,
		INNER_NODE = 0
	};

	int typedIndex;
	int index;
	union
	{
		NodeColor color;
		int height;
	};
	LongType excess;
	LongType capDelta;
	Arc* outArcBegin;
	Arc* outArcEnd;
	Arc* firstOutArc;
	Arc* curArc;
	Arc* parentArc;
	Node* nextInBucket;
	Node* prevInBucket;

	Arc* svArc;
	LongType svCap;
	LongType svSlope;

	Arc* vtArc;
	LongType vtCap;
	LongType vtSlope;
};

struct Arc
{
	LongType flow;
	LongType baseFlow;
	LongType cap;
	Arc* revArc;
	Arc* nextOutArc;

	Node* GetHeadNode();
	Node* GetTailNode();
	
	Node* headNode;

#ifdef VERBOSE
	int index;
#endif
};

struct Bucket
{
	Node* firstNode;
};

struct GraphState
{
	int firstInnerIndex;
	int firstSinkIndex;

	LongType sourceCap;
	LongType sourceSlope;

	LongType sinkCap;
	LongType sinkSlope;

	LongType lambda;
};

typedef std::vector<Node::MinCutType> MinCut;

class Graph
{
public:
	Graph();

	std::vector<Node> nodes;
	std::vector<Arc> arcs;

	GraphState state;
	std::vector<Node*> typedNodeList;

	void SetMultiplier(int multiplier);
	void Reserve(int m, int n);
	void AddAuxArcs(Node* source, Node* sink);
	void PrepareArcs();
	void PrepareNodes(Node* source, Node* sink);

	int GetNodeCount();
	bool IsSourceNode(Node* node);
	bool IsSinkNode(Node* node);
	bool IsInnerNode(Node* node);

	void MakeSource(Node* node);
	void MakeSink(Node* node);
	void RearrangeNodes(MinCut& cut, Node::MinCutType sourceType, int* newFirstinnerIndex, int* newFirstSinkIndex);

	Node* AddNode();
	Arc* AddArcFromSource(Node* tailNode, Node* headNode, int cap, int slope);
	Arc* AddArcToSink(Node* tailNode, Node* headNode, int cap, int slope);
	Arc* AddArc(Node* tailNode, Node* headNode, int cap);

	Node* GetSourceNode();
	Node* GetSinkNode();

	void SaveState(GraphState& state);
	void RestoreState(GraphState& state);

	int ChooseMultiplier(LongType lambda1, LongType lambda2);
	void InitCapDelta(LongType lambda);
	void InitPreflowPush(LongType lambda);
	void IncreaseLambda(LongType lambda);

	LongType GetLowerLambdaBound();
	LongType GetUpperLambdaBound();
	bool FindMaxPreflow(int maxWorkAmount);
	void FindMaxPreflow();
	void Reverse(Graph& revGraph);
	void ConvertPreflowToFlow();
	
	void BackupFlows();
	void FindMinCut(MinCut& minCut);

	void ClearStats();
	void Dump();
	void PrintStats(char* prefix);
	void AccumulateStats(Graph* graph);
	void CheckPreflow();
	void CheckFlow();
	void CheckMaxFlow();
	void InitArcLists();

private:
	void InitCaps(LongType lambda);
	LongType GetArcBaseCap(Arc* arc);
	void DFS(Node* node, std::vector<Node*>& sortedNodes);
	Node* PopNodeForDischarge();
	void DischargeNode(Node* node);
	void PushPreflow(Arc* arc, LongType amount);
	bool RelabelNode(Node* node);
	void ClearBuckets();
	void InitHeights();
	void InitCurArcs();
	void AddToActiveBucket(Node* node);
	void AddToInactiveBucket(Node* node);
	void RemoveFromInactiveBucket(Node* node);
	bool IsArcAlive(Arc* arc);
	void AddAliveArc(Arc* arc);
	void ClearAliveArcs(Node* node);
	void UnmakeSource(Node* node);
	void UnmakeSink(Node* node);
	void Gap(int height);
	void GlobalUpdate();

	std::vector<Bucket> activeBuckets;
	std::vector<Bucket> inactiveBuckets;
	int maxHeight;

	LongType multiplier;
	LongType workSinceLastUpdate;
	LongType workCounter;
	LongType globalUpdateBarrier;
	LongType pushCounter;
	LongType relabelCounter;
	LongType gapCounter;
	LongType globalUpdateCounter;
	LongType cleanupFetchCounter;
	LongType relabelFetchCounter;
};


inline Node* Arc::GetHeadNode()
{
	return headNode;
}

inline Node* Arc::GetTailNode()
{
	return revArc->headNode;
}

}

#endif
