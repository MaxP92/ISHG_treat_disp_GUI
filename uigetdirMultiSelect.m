function [files] = uigetdirMultiSelect( start_path, dialog_title)
% 2018-08 edited by Maxime PINSARD

if nargin == 0 || strcmp(start_path, '') || sum(start_path == 0) % Allow a null argument.
    start_path = pwd;
end

import com.mathworks.mwswing.MJFileChooserPerPlatform;
jchooser = javaObjectEDT('com.mathworks.mwswing.MJFileChooserPerPlatform', start_path);
jchooser.setFileSelectionMode(javax.swing.JFileChooser.DIRECTORIES_ONLY);
jchooser.setMultiSelectionEnabled(true);

if nargin > 1
    jchooser.setDialogTitle(dialog_title);
end

jchooser.showOpenDialog([]);

if jchooser.getState() == javax.swing.JFileChooser.APPROVE_OPTION
    jFiles = jchooser.getSelectedFiles();
    files = arrayfun(@(x) char(x.getPath()), jFiles, 'UniformOutput', false);
elseif jchooser.getState() == javax.swing.JFileChooser.CANCEL_OPTION
    files = [];
else
    error('Error occurred while picking file');
end