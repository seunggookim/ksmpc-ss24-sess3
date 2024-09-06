function [contour] = triouterboundary(surf)
% finds outer boundary of a given 3D trisurface
%
% USAGE
% [contour] = triouterboundary(surf)
%
% INPUT
% surf   (1x1)  triangular mesh structure with .vertices and .faces fields
%
% OUTPUT
% contour (1x1)  structure of which .contour(c).group(g).xyz contains
%                x,y,z-coordinates of g-th contour of c-th level
%
% REF: alecjacobson.com/weblog/?p=1530
%
% (cc) 2020, sgKIM. solleo@gmail.com

% Outerboundary edges are used ONLY once (alecjacobson.com/weblog/?p=1530):
F = surf.faces;
E = sort([F(:,1) F(:,2); F(:,2) F(:,3); F(:,3) F(:,1)]')';
[u,m,n] = unique(E,'rows');
counts = accumarray(n(:), 1);
edges = u(counts==1,:);  % outerboundary edges

% Now connect edges:
isusededge = false(size(edges,1),1);  % is this edge used or not?
g = 1;  % group index
contour = [];
for i = 1:size(edges,1)
  % starting from the i-th edge: e_i = (v_i1, v_i2)
  if ~isusededge(i)
    isusededge(i) = true; % now the i-th edge is used.
    
    verts_thisgroup = edges(i,:);
    % put the vertices of the i-th edge to this group
    
    for v = 1:2 % consider v_i1 and v_i2 separately
      thisv = edges(i,v); % this vertex!
      
      newedges_thisv = find(~~sum(edges == thisv,2) & ~isusededge);
      % this is equivalent to:
      %    setdiff(Edges_including_this_vertex, Edges_used)
      % but much faster!
      
      while ~isempty(newedges_thisv)
        % looping until there's no further connection
        
        newv = edges(newedges_thisv(1),:);
        % just pick one edge from those new edges with "this vertex"
        
        newv(newv == thisv) = [];
        % find the new vertex
        
        if ~isempty(newv)
          % if you really have a new vertex
          
          if v == 1
            verts_thisgroup = [newv verts_thisgroup];
          else
            verts_thisgroup = [verts_thisgroup newv];
          end
          % THIS is the heart of this algorithm. We want to know whether
          % we should put this before or after the current list of
          % vertices to form a one curve.
          
          thisv = newv;
          % now let's look for further connections from the newly added
          % vertex by updating "thisv" with newv
          
          isusededge(newedges_thisv(1)) = true;
          % consider this edge used.
        else
          % this is possible if the edge is looping back (perhaps?)
          
          isusededge(newedges_thisv(1)) = true;
          % just consider this is used and move on.
        end
        
        newedges_thisv = find(~~sum(edges == thisv,2) & ~isusededge);
        % now update new edges with the updated isusededge
        
      end
    end
    contour.group(g).xyz = surf.vertices(verts_thisgroup,:);
    g = g+1;
  end
end
end
