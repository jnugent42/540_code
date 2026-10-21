This is the folder where we are building the automated laser ball characterisation sequence. 

There are three steps to this project:
Step 1: take data cleanly
Step 2: analyse this data accurately
Step 3: combine the analysed data smartly

Step 1: take data cleanly

For each run, before we take it we will want to log what it is including the following information:

Run number
PMT form and theta position for channel 1
PMT form and theta position for channel 2
PMT form and theta position for channel 3
PMT form and theta position for channel 4
is the trigger internal or external 
whether the laser ball is in place
what power laser is being sent 
what is the phi position of the LB
number of samples that we want 
the number of Giga Sample/s
the window width 

To take data run 
python /home/a/usbtmc2/takeDataMain.py 
and follow CAREFULLY the instructions. 

The data will be stored on the lx machines here:
/vols/hyperk/users/your_short_code/WCTE540/LB_Characterisation/data
 make sure you have this folder as well as the 
/vols/hyperk/users/your_short_code/WCTE540/LB_Characterisation/Codes 
  where you'll have copied all of the contents of the /home/a/WCTE540/LB_Characterisation/Codes folder from the 540 lab computer




Multiple things can go bad with this sequence, below are some checks that can be useful to perfom if something looks off:

1) if you do not see anything in the TTS root file: might be that the window where we expect it is wrong, we have set it to 540-650ns in the code TTS_executable.cc which works well when the trigger is set on the ascending slope of the laser syncout with the laser ball. You might have to adjust it if you have a different configuration. try running it once with 0 to 1e3 (if you take 1ms long windows) to pin the position of the TTS peak down.  

2) It can often happen that the PC runs out of space either before a run is complete or before the files are hadded together. If you see that this is the case, make sure you CTRL+C as soon as possible. You can then (after having made space) do python /home/a/usbtmc2/takeDataMain.py again with the same input except for the number of windows which should be put to 0. That way, the code skips the data taking and only does the hadding and combining of existing files.
