function h = rosePatch(cfg)
%     cfg = [];
%     cfg.theta      = PLVang(:,j,3,:);
%     cfg.mode       = 'hist','value'
%     cfg.n_bins     = 90;
%     cfg.rlims      = 700;
%     cfg.axis       = ca;
%     cfg.overlap    = 0;
%     cfg.patch_spec = {'r','facealpha',.3,'facecolor',colors(1,:),'edgecolor',colors(1,:),'linewidth',2};



field = 'rlims'; value = [0 inf];
if ~isfield(cfg,field), cfg.(field) = value; end
field = 'overlap'; value = 0;
if ~isfield(cfg,field), cfg.(field) = value; end

% CREATING AXIS
if isfield(cfg,'axis')
   ca = cfg.axis;
else
   ca = gca;
end

% CREATING AXIS
if isfield(cfg,'axis')
   ca = cfg.axis;
else
   ca = gca;
end

% CALCULATING ANGLES
if isfield(cfg,'complex')
    theta = angle(cfg.complex(:));
    ampl  = abs(cfg.complex(:));
else
    theta = cfg.theta(:);
    ampl  = ones(size(theta));
end
theta = rem(rem(theta,2*pi)+2*pi,2*pi); % Make sure 0 <= theta <= 2*pi

% CALCULATING BINS
if isfield(cfg,'bins')
    bins = sort(rem(cfg.bins(:)',2*pi));
else
    n_bins = cfg.n_bins;
    bins   = (0:n_bins-1)*2*pi/n_bins + pi/n_bins;
end

patch_spec = cfg.patch_spec;

% Determine bin edges and get histogram
edges       = sort(rem([(bins(2:end)+bins(1:end-1))/2 (bins(end)+bins(1)+2*pi)/2],2*pi));
edges       = [edges edges(1)+2*pi];
nn          = histc(rem(theta+2*pi-edges(1),2*pi),edges-edges(1));
nn(end-1)   = nn(end-1)+nn(end);
nn(end)     = [];

% Form radius values for histogram triangle
if min(size(nn))==1, % Vector
  nn = nn(:); 
end
[m,n]       = size(nn);
mm          = 4*m;
r           = zeros(mm,n);
r(2:4:mm,:) = nn;
r(3:4:mm,:) = nn;

% Form theta values for histogram triangle from triangle centers (xx)
zz        = edges;
t         = zeros(mm,1);
t(2:4:mm) = zz(1:m);
t(3:4:mm) = zz(2:m+1);
% Create cirular axis
if numel(cfg.rlims) == 1
    rlims = cfg.rlims;
else
    rlims = min([max([r(:);cfg.rlims(1)]),cfg.rlims(2)]); % setting radius limits
end
if ~cfg.overlap
    polar(ca,0,rlims);
end
% plot the data as patches
[a,b] = pol2cart(t,r);             % convert histogram line to polar coordinates
A     = reshape(a,4,numel(bins));  % reshape 4*N-element x vector into N columns
B     = reshape(b,4,numel(bins));  % reshape 4*N-element y vector into N columns
h     = patch(A,B,patch_spec{:});  % plot N patches based on the columns of A and B

