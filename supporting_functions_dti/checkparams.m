function opts = checkparams(opts_in,lala)
opts = [];
if ~exist('opts_in','var')  || isempty(opts_in)

else
    if mod(size(opts_in,2),2)
        error('option names/values should come in pairs');
    end
    if any(~cellfun(@ischar,opts_in(1:2:end)))
        error('option names should be strings');
    end
    if any(~cellfun(@(x) isstrprop(x(1),'alpha'),opts_in(1:2:end)))
        error('option names should begin with an alphabetic character');
    end
    opt_names = opts_in(1:2:end);
    opt_values = opts_in(2:2:end);
    for i = 1:size(opt_names,2)
        if ischar(opt_values{i})
            opt_values{i} = lower(opt_values{i});
        end
        opts.(lower(opt_names{i})) =  opt_values{i};
    end
end
opts_in = opts;
opts = [];
opts.name = '';
if exist('lala','var')  && ~isempty(lala)
    if mod(size(lala,2),4)
        error('option names/values should come in pairs');
    end
    if any(~cellfun(@ischar,lala(1:4:end)))
        error('option names should be strings');
    end
    if any(~cellfun(@(x) isstrprop(x(1),'alpha'),lala(1:4:end)))
        error('option names should begin with an alphabetic character');
    end
    opt_name  = lala(1:4:end);
    opt_type  = lala(2:4:end);
    opt_range = lala(3:4:end);
    opt_default = lala(4:4:end);
    for i = 1:size(opt_name,2)
        name = opt_name{i};
        type = opt_type{i};
        range = opt_range{i};
        default = opt_default{i};
        if isfield(opts_in,lower(name))
            value = opts_in.(lower(name));
            switch type
                case 'integer'
                    if ischar(value) || ~isscalar(value) || value-round(value)~=0
                        error(['option ' name ' must be an integer']);
                    end
                    if ~isempty(range)
                        if value < range(1) || value > range(2)
                            error(['option ' name ' must be an integer in [' num2str(range(1)) ', ' num2str(range(2)) ']']);
                        end
                    end
                case 'float'
                    disp('float')
                    if ischar(value) || ~isscalar(value)
                        error(['option ' name ' must be a float']);
                    end
                    if ~isempty(range)
                        if value < range(1) || value > range(2)
                            error(['option ' name ' must be a float in [' num2str(range(1)) ', ' num2str(range(2)) ']']);
                        end
                    end
                case 'string'
                    disp('string')
                    if ~ischar(value)
                        error(['option ' name ' must be a string']);
                    end
                    if ~isempty(range)
                        if isempty(intersect(range,value))
                            err = ['option ' name ' must be a string in {'];
                            for j = 1:size(range,2)
                                err = [err '''' range{j} ''', '];
                            end
                            err(end-1) = '}';
                            error(err);
                        end
                    end
                otherwise
                   return;
            end
            opts.(name) = value;
        else
            opts.(name) = default;
        end
    end
end
end
