function [vnew, fnew] = stlDelVerts2(v, f, list)
%STLDELVERT removes a list of vertices from STL files
%V is the Nx3 array of vertices
%F is the Mx3 array of faces
%LIST are the vertices (rows) to delete, where length(LIST) < N for index
%of vertex
%VNEW is the new array of vertices
%FNEW is the new array of faces

% find (on the global set) the position (rows) of the vertices to be deleted
% [~,vdel] = ismember(list,v','rows');

% delete vertices and get new tags
vnew = v;
tags = 1:length(v);
vnew(list,:) = [];
tags(list) = [];

% delete faces
fnew = f;
[fdel,~] = find(ismember(f,list)); % find the position (rows) of the faces to delete
fnew(fdel,:) = [];

% rename faces, as some of the vertices have been deleted
flist = reshape(fnew,numel(fnew),1);
[~,ind] = ismember(flist,tags);
fnew = reshape(ind,numel(flist)/3,3);