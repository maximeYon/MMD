These functions allows for interaction with the ElastiX program for
image registration. 

To get ElastiX running, follow these instructions

1) Download ElastiX from http://elastix.isi.uu.nl

2) Unzip and rename the downloaded folder as 'elastix'

3) Move the program folder into /usr/local/elastix, e.g. by
    'sudo mv ~/Downloads/elastix /usr/local'

4) Typing 'ls /usr/local/elatix' should now show to folder, bin and lib
 
5) To run, you need to be able to start elastix from your terminal by 
   typing 'elastix'. Do this by editing ~/.bashrc, and add the following
   lines 

   export PATH=/usr/local/elastix/bin:/usr/local/elastix/parameters:$PATH
   export DYLD_LIBRARY_PATH=/usr/local/elastix/lib:$DYLD_LIBRARY_PATH
