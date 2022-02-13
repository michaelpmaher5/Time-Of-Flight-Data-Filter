function [ M_edges,E1,E2,E3,E4 ] = edges( S1,S2,S3,S4,Z,nbrOfFrames )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

M_edges = zeros(120,160,nbrOfFrames);
E1=zeros(3,120,nbrOfFrames);
E2=zeros(3,160,nbrOfFrames);
E3=zeros(3,120,nbrOfFrames);
E4=zeros(3,160,nbrOfFrames);


for k=1:nbrOfFrames
    
    % Edge 1
    a1 = (S2(1,1,k)-S1(1,1,k)) / (S2(2,1,k)-S1(2,1,k));
    b1 = S1(1,1,k) - a1 * S1(2,1,k);
    last_valid_i = S1(2,1,k);
    last_valid_j = S1(1,1,k);
    Holes = [];
    
    
    for i=S1(2,1,k):S2(2,1,k)
        
        j = round(a1*i+b1);
        right_limit = j+3;
        
        while Z(i,j,k)==0 && j<right_limit
            j=j+1;
        end
        
        if j<right_limit
            % No hole
            E1(:,i,k) = [j; i; Z(i,j,k)];
            
            if ~isempty(Holes)
                % Fixing the previous holes
                a = (last_valid_j-j) / (last_valid_i-i);
                b = j - a * i;
                for i_hole=Holes
                    j_hole = round(a*i_hole+b);
                    if Z(i_hole,j_hole,k)~=0
                        E1(:,i_hole,k) = [j_hole; i_hole; Z(i_hole,j_hole,k)];
                    else
                        E1(:,i_hole,k) = [j_hole; i_hole; (Z(i,j,k)*norm([j_hole-last_valid_j i_hole-last_valid_i])+Z(last_valid_i,last_valid_j,k)*norm([j-j_hole i-i_hole]))/norm([j-last_valid_j i-last_valid_i])];
                    end
                end
                Holes = [];
            end
            
            last_valid_i = i;
            last_valid_j = j;
            
        else
            % Hole
            Holes = [Holes i];
        end
        
    end
    
    
    % Edge 2
    a2 = (S3(2,1,k)-S2(2,1,k)) / (S3(1,1,k)-S2(1,1,k));
    b2 = S2(2,1,k) - a2 * S2(1,1,k);
    last_valid_i = S2(2,1,k);
    last_valid_j = S2(1,1,k);
    Holes = [];
    
    
    for j=S2(1,1,k):S3(1,1,k)
        
        i = round(a2*j+b2);
        down_limit = i-3;
        
        while Z(i,j,k)==0 && i>down_limit
            i=i-1;
        end
        
        if i>down_limit
            % No hole
            E2(:,j,k) = [j; i; Z(i,j,k)];
            
            if ~isempty(Holes)
                % Fixing the previous holes
                a = (last_valid_i-i) / (last_valid_j-j);
                b = i - a * j;
                for j_hole=Holes
                    i_hole = round(a*j_hole+b);
                    if Z(i_hole,j_hole,k)~=0
                        E2(:,j_hole,k) = [j_hole; i_hole; Z(i_hole,j_hole,k)];
                    else
                        E2(:,j_hole,k) = [j_hole; i_hole; (Z(i,j,k)*norm([j_hole-last_valid_j i_hole-last_valid_i])+Z(last_valid_i,last_valid_j,k)*norm([j-j_hole i-i_hole]))/norm([j-last_valid_j i-last_valid_i])];
                    end
                end
                Holes = [];
            end
            
            last_valid_i = i;
            last_valid_j = j;
            
        else
            % Hole
            Holes = [Holes j];
        end
        
    end
    
    
    % Edge 3
    a3 = (S4(1,1,k)-S3(1,1,k)) / (S4(2,1,k)-S3(2,1,k));
    b3 = S3(1,1,k) - a3 * S3(2,1,k);
    last_valid_i = S3(2,1,k);
    last_valid_j = S3(1,1,k);
    Holes = [];
    
    
    for i=S3(2,1,k):-1:S4(2,1,k)
        
        j = round(a3*i+b3);
        left_limit = j-3;
        
        while Z(i,j,k)==0 && j>left_limit
            j=j-1;
        end
        
        if j>left_limit
            % No hole
            E3(:,i,k) = [j; i; Z(i,j,k)];
            
            if ~isempty(Holes)
                % Fixing the previous holes
                a = (last_valid_j-j) / (last_valid_i-i);
                b = j - a * i;
                for i_hole=Holes
                    j_hole = round(a*i_hole+b);
                    if Z(i_hole,j_hole,k)~=0
                        E3(:,i_hole,k) = [j_hole; i_hole; Z(i_hole,j_hole,k)];
                    else
                        E3(:,i_hole,k) = [j_hole; i_hole; (Z(i,j,k)*norm([j_hole-last_valid_j i_hole-last_valid_i])+Z(last_valid_i,last_valid_j,k)*norm([j-j_hole i-i_hole]))/norm([j-last_valid_j i-last_valid_i])];
                    end
                end
                Holes = [];
            end
            
            last_valid_i = i;
            last_valid_j = j;
            
        else
            % Hole
            Holes = [Holes i];
        end
        
    end
    
    
    % Edge 4
    a4 = (S1(2,1,k)-S4(2,1,k)) / (S1(1,1,k)-S4(1,1,k));
    b4 = S4(2,1,k) - a4 * S4(1,1,k);
    last_valid_i = S4(2,1,k);
    last_valid_j = S4(1,1,k);
    Holes = [];
    
    
    for j=S4(1,1,k):-1:S1(1,1,k)
        
        i = round(a4*j+b4);
        up_limit = i+3;
        
        while Z(i,j,k)==0 && i<up_limit
            i=i+1;
        end
        
        if i<up_limit
            % No hole
            E4(:,j,k) = [j; i; Z(i,j,k)];
            
            if ~isempty(Holes)
                % Fixing the previous holes
                a = (last_valid_i-i) / (last_valid_j-j);
                b = i - a * j;
                for j_hole=Holes
                    i_hole = round(a*j_hole+b);
                    if Z(i_hole,j_hole,k)~=0
                        E4(:,j_hole,k) = [j_hole; i_hole; Z(i_hole,j_hole,k)];
                    else
                        E4(:,j_hole,k) = [j_hole; i_hole; (Z(i,j,k)*norm([j_hole-last_valid_j i_hole-last_valid_i])+Z(last_valid_i,last_valid_j,k)*norm([j-j_hole i-i_hole]))/norm([j-last_valid_j i-last_valid_i])];
                    end
                end
                Holes = [];
            end
            
            last_valid_i = i;
            last_valid_j = j;
            
        else
            % Hole
            Holes = [Holes j];
        end
        
    end
    
    E1(1,:,k) = round(avg_7(E1(1,:,k)));
    E2(2,:,k) = round(avg_7(E2(2,:,k)));
    E3(1,:,k) = round(avg_7(E3(1,:,k)));
    E4(2,:,k) = round(avg_7(E4(2,:,k)));
    
    for p=find(E1(1,:,k))
        M_edges(E1(2,p,k),E1(1,p,k),k) = E1(3,p,k);
    end
    for p=find(E2(1,:,k))
        M_edges(E2(2,p,k),E2(1,p,k),k) = E2(3,p,k);
    end
    for p=find(E3(1,:,k))
        M_edges(E3(2,p,k),E3(1,p,k),k) = E3(3,p,k);
    end
    for p=find(E4(1,:,k))
        M_edges(E4(2,p,k),E4(1,p,k),k) = E4(3,p,k);
    end
    
    E1(1,:,k) = avg_7(E1(1,:,k));
    E2(2,:,k) = avg_7(E2(2,:,k));
    E3(1,:,k) = avg_7(E3(1,:,k));
    E4(2,:,k) = avg_7(E4(2,:,k));

end

end

