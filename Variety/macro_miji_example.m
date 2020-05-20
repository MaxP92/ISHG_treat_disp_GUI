%
% Working example run macro in Matlab using Miji 
% 2017.5.16
% MatLab R2016a and Fiji 1.51n 64bits

javaaddpath 'C:\Program Files\MATLAB\R2016a\java\jar\mij.jar'
javaaddpath 'C:\Program Files\MATLAB\R2016a\java\jar\ij-1.51n.jar'

Miji(false);

IJ=ij.IJ(); 
macro_path = 'C:\Users\pc\Documents\These\Biomedical\CD SHG\codes CD\Macro ImageJ CD\cd_Z.ijm';
IJ.runMacroFile(java.lang.String(macro_path),java.lang.String()); 

a=magic(25);
MIJ.createImage('result', a, true);

MIJ.run('StackReg', 'transformation=Translation');



