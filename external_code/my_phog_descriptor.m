function p = my_phog_descriptor(bh,bv,L,bin)
% anna_PHOGDESCRIPTOR Computes Pyramid Histogram of Oriented Gradient over a ROI.
%               
%IN:
%	bh - matrix of bin histogram values
%	bv - matrix of gradient values 
%   L - number of pyramid levels
%   bin - number of bins
%
%OUT:
%	p - pyramid histogram of oriented gradients (phog descriptor)
if(isempty(bh) && isempty(bv))
  error('no edges!');
end

if(L==0)
    size_feat = 1*bin;
elseif(L==1)
    size_feat = bin + 4*bin;
elseif(L==2)    
    size_feat = bin + 4*bin + 16*bin;
elseif(L==3)    
    size_feat = bin + 4*bin + 16*bin + 64* bin;
else
    error('not prepared for that');
end
    
%level 0

duh = cell(bin,1);
%duh2 = zeros(size(bv,1),size(bv,2), bin);
p = zeros(size_feat,1);

counter = 0;
for b=1:bin
        
    duh{b} = bv.*(bh==b);    
    p(counter + 1) = sum(sum(duh{b}));
    
    %duh2(:,:,b) = bv.*(bh==b);
    %p(counter + 1) = sum(sum(duh2(:,:,b)));    
    %assert(all(p2 == p));
    %assert(all(all(duh{b} == duh2(:,:,b))))
    counter = counter + 1;
end        

cella = 1;

for l=1:L 
    x = fix(size(bh,2)/(2^l));
    y = fix(size(bh,1)/(2^l));
    xx=0;
    yy=0;
    
    counterx = 0;
    while xx+x<=size(bh,2)        
       countery=0;
       while yy +y <=size(bh,1) 
            
            %bh_cella = bh(yy+1:yy+y,xx+1:xx+x);
            %bv_cella = bv(yy+1:yy+y,xx+1:xx+x);                        
                        
            for b=1:bin
                %p(counter+1) = sum(bv_cella(bh_cella==b));
                
                p(counter+1) = sum(sum(duh{b}(yy+1:yy+y,xx+1:xx+x)));
                %p(counter+1) = sum(sum(duh2(yy+1:yy+y,xx+1:xx+x, b)));
                counter = counter + 1;
            end 
            yy = yy+y;
            countery = countery+1;
            if(countery == (2^l))
                break;
            end
        end        
        cella = cella+1;
        yy = 0;
        xx = xx+x;
        
        counterx = counterx + 1;
        if(counterx == (2^l))
            break;
        end
    end
end
% if sum(p)~=0
%     p = p/sum(p);
% end

