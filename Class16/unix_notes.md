# Basic Unix 

Some important file system commands include 

pwd: Print working directory
ls: List files and folders
cd: Change directory
mkdir: Make new directory 
rm: Delete files and directories
nano: A very basic text editor that is always available
less: To view/read text files page by page(pager program)

My AWS instance: 
ssh -i ~/Downloads/bimm143_kaliyah.pem ubuntu@ec2-44-246-247-67.us-west-2.compute.amazonaws.com

To copy from my AWS instance:
scp -i ~/Downloads/bimm143_kaliyah.pem ubuntu@ec2-44-246-247-67.us-west-2.compute.amazonaws.com:~/work/results.tsv .

open New Terminal (commands)
- whoami
- cd ~/Desktop/class16
- pwd
- scp -i ~/Downloads/bimm143_kaliyah.pem ubuntu@ec2-44-246-247-67.us-west-2.compute.amazonaws.com:~/work/results.tsv .
- ls