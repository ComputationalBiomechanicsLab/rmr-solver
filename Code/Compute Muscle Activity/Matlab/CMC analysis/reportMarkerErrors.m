function reportMarkerErrors(dirName)

trials ={ 'abd01'; 'abd02'; 'abd03';...
            'abd21'; 'abd22'; 'abd23';...
            'flx01'; 'flx02'; 'flx03';...
            'flx21'; 'flx22'; 'flx23';...
            'shrug01'; 'shrug02'; 'shrug03';...
            'shrug21'; 'shrug22'; 'shrug23' };

if(nargin < 1)
	dirName = 'output';
end

for i=1:18
    file=[dirName '/' trials{i} '_ik_marker_errors.sto'];
    [time,squared,rms,max_e]=textread(file,'%f %f %f %f','headerlines',7);
    rms_mean(i)=mean(rms);
    max_e_max(i)=max(max_e);
end

fprintf('%s\t', trials{:});
fprintf('\n');
fprintf('%1.4f\t', rms_mean);
fprintf('\n');
fprintf('%1.4f\t', max_e_max);
fprintf('\n');
fprintf('Mean of all RMSE: %g', mean(rms_mean));
fprintf('\n');
fprintf('Max of max error: %g', max(max_e_max));
fprintf('\n');
