%demo_distmesh

%% 3D tet

%% Example: Unit ball
 dim=3;
%  d=inline('(sum(p.^2,2))-1','p');
%  [p,t]=distmeshnd(d,@huniform,0.2,[-ones(1,dim);ones(1,dim)],[]);

 fd=inline('sqrt(sum(p.^2,2))-1','p');
[p,t]=distmeshnd(fd,@huniform,0.15,[-1,-1,-1;1,1,1],[]);
%% Example: (Uniform Mesh on Unit Sphere)
fd=@(p) dsphere(p,0,0,0,1);
[p,t]=distmeshsurface(fd,@huniform,0.2,1.1*[-1,-1,-1;1,1,1]);
    
%% Example: (Graded Mesh on Unit Sphere)
fd=@(p) dsphere(p,0,0,0,1);
fh=@(p) 0.05+0.5*dsphere(p,0,0,1,0);
[p,t]=distmeshsurface(fd,fh,0.15,1.1*[-1,-1,-1;1,1,1]);
    
%% Example: (Uniform Mesh on Torus)
fd=@(p) (sum(p.^2,2)+.8^2-.2^2).^2-4*.8^2*(p(:,1).^2+p(:,2).^2);
[p,t]=distmeshsurface(fd,@huniform,0.1,[-1.1,-1.1,-.25;1.1,1.1,.25]);

%% Example: (Uniform Mesh on Ellipsoid)
fd=@(p) p(:,1).^2/4+p(:,2).^2/1+p(:,3).^2/1.5^2-1;
[p,t]=distmeshsurface(fd,@huniform,0.2,[-2.1,-1.1,-1.6; 2.1,1.1,1.6]);
  


%% 2D


%% Example: (Uniform Mesh on Unit Circle)
    fd=@(p) sqrt(sum(p.^2,2))-1;
    [p,t]=distmesh2d(fd,@huniform,0.2,[-1,-1;1,1],[]);
%% Example: (Rectangle with circular hole, refined at circle boundary)
    fd=@(p) ddiff(drectangle(p,-1,1,-1,1),dcircle(p,0,0,0.5));
    fh=@(p) 0.05+0.3*dcircle(p,0,0,0.5);
    [p,t]=distmesh2d(fd,fh,0.05,[-1,-1;1,1],[-1,-1;-1,1;1,-1;1,1]);
%% Example: (Polygon)
pv=[-0.4 -0.5;0.4 -0.2;0.4 -0.7;1.5 -0.4;0.9 0.1;
    1.6 0.8;0.5 0.5;0.2 1;0.1 0.4;-0.7 0.7;-0.4 -0.5];
[p,t]=distmesh2d(@dpoly,@huniform,0.1,[-1,-1; 2,1],pv,pv);

%% Example: (Ellipse)
fd=@(p) p(:,1).^2/2^2+p(:,2).^2/1^2-1;
[p,t]=distmesh2d(fd,@huniform,0.2,[-2,-1;2,1],[]);
    
%% Example: (Square, with size function point and line sources)
fd=@(p) drectangle(p,0,1,0,1);
fh=@(p) min(min(0.01+0.3*abs(dcircle(p,0,0,0)), ...
             0.025+0.3*abs(dpoly(p,[0.3,0.7; 0.7,0.5]))),0.15);
[p,t]=distmesh2d(fd,fh,0.01,[0,0;1,1],[0,0;1,0;0,1;1,1]);




%% Example: (NACA0012 airfoil)
hlead=0.01; htrail=0.04; hmax=2; circx=2; circr=4;
a=.12/.2*[0.2969,-0.1260,-0.3516,0.2843,-0.1036];

fd=@(p) ddiff(dcircle(p,circx,0,circr),(abs(p(:,2))-polyval([a(5:-1:2),0],p(:,1))).^2-a(1)^2*p(:,1));
fh=@(p) min(min(hlead+0.3*dcircle(p,0,0,0),htrail+0.3*dcircle(p,1,0,0)),hmax);

fixx=1-htrail*cumsum(1.3.^(0:4)');
fixy=a(1)*sqrt(fixx)+polyval([a(5:-1:2),0],fixx);
fix=[[circx+[-1,1,0,0]*circr; 0,0,circr*[-1,1]]'; 0,0; 1,0; fixx,fixy; fixx,-fixy];
box=[circx-circr,-circr; circx+circr,circr];
h0=min([hlead,htrail,hmax]);

[p,t]=distmesh2d(fd,fh,h0,box,fix);