This folder collects MATLAB scripts commonly used for processing raw data.

In `OpenSim4.3` you can find scripts that use OpenSim functionalities specifically. They described here, to view all of them on the same page

When adding a new file, please edit the following list with the file name and a short description of how to use it (and why).

- `EMGcrunch_example.m` : example of how to process EMG signals of specific (user defined muscles). Requires MVC experiments and generic motions to be analyzed (in format `MVC*.mat` and `motion*.mat`). The EMG data are filtered, rectified, and normalized.

- `OpenSim4.3/c3dExport_markersOnly.m` : extracts marker trajectories (`.trc`) from motion capture results (`.c3d`). When run, prompts the user to navigate to the desired `.c3d` file, processes it and saves it with a user-defined name in the same folder.

- `EMG_analysis.m`: processes EMG signals of specific (user defined muscles). Requires MVC experiments and generic motions to be analyzed (in format `MVC*.mat` and `motion*.mat`). The EMG data are filtered, rectified, and normalized.
		    Uses: `butterworth_filter.m`, `extract_data.m`, `extract_mvc.m`, and `normalize_emg.m`

- `butterworth_filter.m`: filters raw data by applying a bandpass filter, rectification and low pass filter to the EMG data. Output is filtered EMG prepped for analysis

- `extract_data.m`: takes `.mat` files with analog data as input and extracts EMG channels. NB: EMG channel should be the first channel per muscle.

- `extract_mvc.m` takes the maximum of the maximum of the filtered data, this results in one value per muscle: the highest measured activation over all trials.

- `normalize_emg.m`: normalizes trial data with the use of the MVC data. Scales activation between 0 and 1.

- `readTRC.m`: reads a `.trc` file (only input) and loads it into the workspace. It returns the x-y-z trajectories for each marker in form of a single 2D array. Timestamps, labels for the markers and units of measurement used can also be retrieved.

- `writeMarkersToTRC.m`: saves a new `.trc` file. Takes as inputs the name of the files, the markers to be saved (3D markers trajectories as a single 2D array), marker labels, rate of acquisition, optionally the frames, the timestamps and the units of measurement used. A GUI will be prompted at the end to ask the directory for the `.trc' file to be saved in.

- `read_motionFile.m`: reads a `.mot` file (only input) and returns a Matlab structure with fields corresponding to column labels, 2D matrix of data, dimensions of the matrix of data.

- `CSV2TRC.m`: opens a `.CSV` file (or other ASCII text file that can be converted to a table) and transforms it to match the `.trc` format. Uses `writeMarkersToTRC.m`.


- ... [add new script here]

