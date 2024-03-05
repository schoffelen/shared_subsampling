function run_script(scriptname, varargin)

if numel(varargin)>0
  for k = 1:numel(varargin)
    eval([varargin{k }{1},'=varargin{k}{2}']);
  end
end
eval(scriptname);
