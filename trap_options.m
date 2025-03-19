function [ncol,m,reltol,itmax,verbose,gmrestol,spurioustol,lusetting] = trap_options(varargin)

ncol     = 200;
m        = [8,8];
reltol   = 1e-8;
itmax    = 10;
verbose  = 0;
gmrestol = min(1e-13,reltol/10);
spurioustol = 1e-3;
lusetting = [0.1 0.001];
if nargin == 1 && ~isnumeric(varargin{1})
    opt = varargin{1};
    if isfield(opt,'ncol')
        ncol = opt.ncol;
    end
    if isfield(opt,'m')
        m = opt.m;
    end
    if isfield(opt,'reltol')
        reltol = opt.reltol;
    end
    if isfield(opt,'itmax')
        itmax = opt.itmax;
    end
    if isfield(opt,'verbose')
        verbose = opt.verbose;
    end
    if isfield(opt,'gmrestol')
        gmrestol = opt.gmrestol;
    end
    if isfield(opt,'spurioustol')
        spurioustol = opt.spurioustol;
    end
    if isfield(opt,'lusetting')
        lusetting = opt.lusetting;
    end
end

if nargin >= 1 && isnumeric(varargin{1})
    if nargin >= 1
        ncol = varargin{1};
    end
    if nargin >= 2
        m = varargin{2};
    end
    if nargin >= 3
        reltol = varargin{3};
    end
    if nargin >= 4
        itmax = varargin{4};
    end
    if nargin >= 5
        verbose = varargin{5};
    end
    if nargin >= 6
        gmrestol = varargin{6};
    end
    if nargin >= 7
        spurioustol = varargin{7};
    end
    if nargin >= 8
        lusetting = varargin{8};
    end
end
end