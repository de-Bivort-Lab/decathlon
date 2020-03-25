function cross_validate_pcs(data,varargin)


plot_bool = true;
test_size = 0.2;
for i=1:numel(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            case 'Plot'
                i=i+1;
                plot_bool = varargin{i};
            case 'TestSize'
                i=i+1;
                test_size = varargin{i};
        end
    end
end