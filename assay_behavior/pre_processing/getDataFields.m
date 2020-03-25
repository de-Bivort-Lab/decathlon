function [data,f] = getDataFields(expmt)

    reset(expmt);
    if isfield(expmt.data,'speed')
        speed = expmt.data.speed.raw();
        log_speed = log(speed(:));
        [~,mus] = arrayfun(@(i) kthresh_distribution(log_speed), 1:20,...
                'UniformOutput', false);
        mus = exp(nanmedian(cat(2,mus{:}),2));
        thresh = mus(1);
        log_thresh = log(thresh);
    elseif isfield(expmt.meta,'velocity')
        speed = expmt.meta.velocity;
        log_speed = log(speed(:));
        [~,mus] = arrayfun(@(i) kthresh_distribution(log_speed), 1:20,...
                'UniformOutput', false);
        mus = exp(nanmedian(cat(2,mus{:}),2));
        thresh = mus(1);
        log_thresh = log(thresh);
    end


    switch expmt.meta.name
        
        case 'Arena'
            
            f = {'circling';'speed'}; 
            data.circling = expmt.meta.handedness.mu;
            data.speed = nanmean(speed);
            data.filter = data.speed > thresh;

        case 'Y-Maze'     
            
            idx = 1:numel(expmt.meta.labels_table.ID);
            expmt.meta.num_traces = numel(expmt.meta.labels_table.ID);
            
            f = {'circling_mu';'right_bias';'hand_clumpiness';...
                'hand_switchiness';'speed';'nTrials'};
            data.circling_mu = expmt.meta.handedness.mu(idx);
            data.right_bias = expmt.data.Turns.rBias(idx);
            %data.hand_clumpiness = expmt.data.Turns.clumpiness(idx);
            data.hand_switchiness = expmt.data.Turns.switchiness(idx);
            data.speed = nanmean(expmt.data.speed.raw(:,idx));
            data.nTrials = expmt.data.Turns.n(idx);
            data.filter = data.nTrials > 10;

            % re-calculate clumpiness
            t_idx = num2cell(expmt.data.Turns.raw()~=0,1);
            t = cumsum(expmt.data.time.raw());
            iti = cellfun(@(i) diff(t(i)), t_idx, 'UniformOutput', false);
            n = expmt.data.Turns.n;
            n(n<0) = 0;
            clumpiness = cellfun(@(i,n) std(i)/(t(end)/(n+1)), iti, num2cell(n));
            data.hand_clumpiness = clumpiness(idx);

        case 'LED Y-maze'

            idx = 1:numel(expmt.meta.labels_table.ID);
            expmt.meta.num_traces = numel(expmt.meta.labels_table.ID);
            
            f = {'circling_mu';'right_bias';'light_bias';'hand_clumpiness';...
                'hand_switchiness';'light_switchiness';'speed';'nTrials'}; 
            data.circling_mu = expmt.meta.handedness.mu(idx);
            data.right_bias = expmt.data.Turns.rBias(idx);
            data.light_bias = expmt.data.LightChoice.pBias(idx);
            data.speed = nanmean(expmt.data.speed.raw(:,idx));
            data.nTrials = expmt.data.Turns.n(idx);
            data.hand_switchiness = expmt.data.Turns.switchiness(idx);
            data.light_switchiness = expmt.data.LightChoice.switchiness(idx);
            data.filter = data.nTrials > 10;

            % re-calculate clumpiness
            t_idx = num2cell(expmt.data.Turns.raw()~=0,1);
            t = cumsum(expmt.data.time.raw());
            iti = cellfun(@(i) diff(t(i)), t_idx, 'UniformOutput', false);
            n = expmt.data.Turns.n;
            n(n<0) = 0;
            clumpiness = cellfun(@(i,n) std(i)/(t(end)/(n+1)), iti, num2cell(n));
            data.hand_clumpiness = clumpiness(idx);

        case 'Slow Phototaxis'

            f = {'circling';'speed';'occupancy';'nTrials'}; 
            data.circling = expmt.meta.handedness.mu;
            data.occupancy = expmt.meta.Light.avg_occ;
            data.speed = nanmean(speed);
            data.nTrials = cellfun(@(t) sum(t>0), expmt.meta.Light.tInc);
            data.filter = data.speed > thresh;

        case 'Optomotor'

            f = {'circling';'speed';'optomotor_index';'nTrials'}; 
            data.circling = expmt.meta.handedness.mu;
            data.optomotor_index = -expmt.meta.Optomotor.index;
            data.speed = nanmean(speed);
            data.nTrials = sum(diff(expmt.data.StimStatus.raw())==1);
            data.filter = data.speed > thresh;
            
        case 'Circadian'
            
            f = {'circling';'speed';'gravitactic_index'}; 
            data.circling = expmt.meta.handedness.mu;
            data.gravitactic_index = expmt.data.area.gravity_index;
            data.speed = nanmean(speed);
            data.filter = data.speed > thresh;

        case 'Olfaction'
            
            f = {'occupancy';'right_bias';'hand_clumpiness';...
                'hand_switchiness';'speed';'nTrials';'preodor_occupancy'};
            data.occupancy = expmt.meta.occupancy;
            data.right_bias = expmt.data.Turns.rBias;
            data.hand_clumpiness = expmt.data.Turns.clumpiness;
            data.hand_switchiness = expmt.data.Turns.switchiness;
            data.speed = nanmean(speed);
            data.nTrials = expmt.data.Turns.n;
            data.preodor_occupancy = expmt.meta.pre_odor_occupancy;
            data.filter = data.speed > thresh;
            
        case 'Temporal Phototaxis'
            
            f = {'circling';'speed';'occupancy';'light_clumpiness';'nTrials'}; 
            data.circling = expmt.meta.handedness.mu;
            data.speed = nanmean(speed);
            data.occupancy = expmt.data.LightStatus.occ;   
            data.nTrials = expmt.data.LightStatus.n;
            data.filter = data.speed > thresh;

            % re-calculate clumpiness
            t = cumsum(expmt.data.time.raw());
            iti = cellfun(@(i) diff([0; t(i)]), expmt.data.LightStatus.trans,...
                'UniformOutput', false);
            n = cellfun(@numel,expmt.data.LightStatus.trans);
            n(n<0) = 0;
            data.light_clumpiness = cellfun(@(i,n) std(i)/(t(end)/(n+1)), iti, num2cell(n));

        otherwise
            errordlg('Experiment name not recognized, no analysis performed');
            
    end

    if isfield_all(expmt,'data.speed.bouts.idx')
        % add new movement bout fields
        f = [f; {'nBouts';'bout_length';'bout_clumpiness'}];

        % calculate num bouts and bout length
        if exist('idx','var')
            idx = expmt.data.speed.bouts.idx(idx);
        else
            idx = expmt.data.speed.bouts.idx;
        end
        idx(cellfun(@isempty,idx)) = {zeros(0,2)};
        data.nBouts = cellfun(@(i) size(i,1), idx);
        data.bout_length = cellfun(@(i) nanmean(diff(i,1,2)),idx);

        % calculate bout length clumpiness
        t = cumsum(expmt.data.time.raw());
        dur = t(end);
        iti = cellfun(@(i) diff([[0; t(i(1:end-1,2))] t(i(:,1))],1,2), idx,...
            'UniformOutput', false);
        data.bout_clumpiness = cellfun(@(i,n) std(i)/(dur/(n+1)), iti, num2cell(data.nBouts));
    end
    
    % standardize dimensions
    fn = fieldnames(data);
    for i = 1:length(fn)
        tmp = data.(fn{i});
        if find(size(tmp)==expmt.meta.num_traces,1)==2
            data.(fn{i}) = data.(fn{i})';
        end
    end      
        
    
end