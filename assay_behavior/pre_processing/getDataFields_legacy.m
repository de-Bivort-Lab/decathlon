function [data,fn] = getDataFields_legacy(expmt)

    if isfield(expmt,'Speed') && isfield(expmt.Speed,'data')
        speed = expmt.Speed.data;
        log_speed = log(speed(:));
        [~,mus] = arrayfun(@(i) kthresh_distribution(log_speed), 1:20,...
                'UniformOutput', false);
        mus = exp(nanmedian(cat(2,mus{:}),2));
        thresh = mus(1);
    elseif isfield(expmt,'velocity')
        speed = expmt.velocity;
        log_speed = log(speed(:));
        [~,mus] = arrayfun(@(i) kthresh_distribution(log_speed), 1:20,...
                'UniformOutput', false);
        mus = exp(nanmedian(cat(2,mus{:}),2));
        thresh = mus(1);
    end

    switch expmt.Name
        
        case 'Arena'
            
            data.circling = expmt.handedness.mu;
            data.speed = nanmean(speed);
            data.filter = data.speed > thresh;

        case 'Y-maze'     
            
            idx = 1:size(expmt.labels_table,1);
            expmt.nTracks = length(idx);
            
            data.circling = expmt.handedness.mu(idx);
            data.right_bias = expmt.Turns.rBias(idx);
            data.nTrials = expmt.Turns.n(idx);
            data.hand_clumpiness = expmt.Turns.clumpiness(idx);
            data.hand_switchiness = expmt.Turns.switchiness(idx);
            data.speed = nanmean(speed(:,idx));
            data.filter = data.nTrials > 10;

            % re-calculate clumpiness
            t_idx = num2cell(~isnan(expmt.Turns.data),1);
            t = cumsum(expmt.Time.data);
            iti = cellfun(@(i) [0;diff(t(i))], t_idx, 'UniformOutput', false);
            n = expmt.Turns.n;
            n(n<0) = 0;
            clumpiness = cellfun(@(i,n) std(i)/(t(end)/(n+1)), iti, num2cell(n));
            data.hand_clumpiness = clumpiness(idx);

        case 'LED Y-maze'

            idx = expmt.labels{1,4}:expmt.labels{1,5};
            expmt.nTracks = length(idx);

            data.circling = expmt.handedness.mu(idx);
            data.right_bias = expmt.Turns.rBias(idx);
            data.light_bias = expmt.LightChoice.pBias(idx);
            data.nTrials = expmt.Turns.n(idx);
            data.hand_clumpiness = expmt.Turns.clumpiness(idx);
            data.hand_switchiness = expmt.Turns.switchiness(idx);
            data.light_switchiness = expmt.LightChoice.switchiness(idx);
            data.speed = nanmean(speed(:,idx));
            data.filter = data.nTrials > 10;

            % re-calculate clumpiness
            t_idx = num2cell(~isnan(expmt.Turns.data),1);
            t = cumsum(expmt.Time.data);
            iti = cellfun(@(i) [0;diff(t(i))], t_idx, 'UniformOutput', false);
            n = expmt.Turns.n;
            n(n<0) = 0;
            clumpiness = cellfun(@(i,n) std(i)/(t(end)/(n+1)), iti, num2cell(n));
            data.hand_clumpiness = clumpiness(idx);

        case 'Slow Phototaxis'

            data.circling = expmt.handedness.mu;
            data.circling_blank = expmt.handedness_Blank.mu;
            data.nTrials = sum(cell2mat(expmt.Light.tInc)>0.005);
            data.occupancy = expmt.Light.avg_occ;
            data.speed = nanmean(speed);
            data.filter = data.speed > thresh;

        case 'Optomotor'

            data.circling = expmt.handedness.mu;
            data.optomotor_index = -expmt.Optomotor.index;
            data.speed = nanmean(speed);
            data.nTrials = expmt.Optomotor.n;
            data.filter = data.speed > thresh;

        case 'Circadian'

            data.circling = expmt.handedness.mu;
            data.speed = nanmean(speed);
            data.filter = data.speed > thresh;
            if isfield(expmt,'Gravity')
                data.gravitactic_index = expmt.Gravity.index;
            end
            
        case 'Temporal Phototaxis'
            
            data.circling = expmt.handedness.mu;
            data.occupancy = expmt.LightStatus.occ;
            data.nTrials = expmt.Lightstatus.n;
            data.speed = nanmean(speed);
            data.filter = data.speed > thresh;

            t = cumsum(expmt.Time.data);
            iti = cellfun(@(i) diff([0; t(i)]), expmt.LightStatus.trans,...
                'UniformOutput', false);
            n = cellfun(@numel,expmt.LightStatus.trans);
            n(n<0) = 0;
            data.light_clumpiness = cellfun(@(i,n) std(i)/(t(end)/(n+1)), iti, num2cell(n));
            
        case 'Olfaction'
            
            data.occupancy = expmt.occupancy;
            data.preodor_occupancy = expmt.preOdorOccupancy;
            data.right_bias = expmt.Turns.rBias;
            data.hand_clumpiness = expmt.Turns.clumpiness;
            data.hand_switchiness = expmt.Turns.switchiness;
            data.nTrials = expmt.Turns.n;
            data.speed = nanmean(speed);
            data.filter = data.speed > thresh;

        otherwise
            errordlg('Experiment name not recognized, no analysis performed');
    end

    if isfield_all(expmt,'Speed.bout_idx')

        % calculate num bouts and bout length
        if exist('idx','var')
            idx = expmt.Speed.bout_idx(idx);
        else
            idx = expmt.Speed.bout_idx;
        end
        idx(cellfun(@isempty,idx)) = {zeros(0,2)};
        data.nBouts = cellfun(@(i) size(i,1), idx);
        data.bout_length = cellfun(@(i) nanmean(diff(i,1,2)),idx);

        % calculate bout length clumpiness
        t = cumsum(expmt.Time.data);
        dur = t(end);
        
        iti = cellfun(@(i) diff([[0; t(i(1:end-1,2))] t(i(:,1))],1,2), idx,...
            'UniformOutput', false);
        data.bout_clumpiness = cellfun(@(i,n) std(i)/(dur/(n+1)), iti, num2cell(data.nBouts));
    end
    
    % standardize dimensions
    fn = fieldnames(data);
    for i = 1:length(fn)
        tmp = data.(fn{i});
        if find(size(tmp)==expmt.nTracks,1)==2
            data.(fn{i}) = data.(fn{i})';
        end
    end
    
    fn(strmatch('filter',fn))=[];
        
    
end