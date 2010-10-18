###
# Align the data with Bowtie and dump the resulting BAM files
# on S3
# Run Macs peak finding on them.

require 'catpaws'

#generic settings
set :aws_access_key,  ENV['AMAZON_ACCESS_KEY']
set :aws_secret_access_key , ENV['AMAZON_SECRET_ACCESS_KEY']
set :ec2_url, ENV['EC2_URL']
set :ssh_options, { :user => "ubuntu", :keys=>[ENV['EC2_KEYFILE']]}
set :key, ENV['EC2_KEY']
set :key_file, ENV['EC2_KEYFILE']
set :ami, 'ami-20794c54'  #EC2 eu-west-1 64bit Lucid
set :instance_type, 'm1.large'
set :s3cfg, ENV['S3CFG'] #location of ubuntu s3cfg file
set :working_dir, '/mnt/work'

#note to self
#ami-52794c26 32bit Lucid
#ami-505c6924 64bit Maverick
#ami-20794c54 64bit Lucid

set :nhosts, 1
set :group_name, 'ns5dastro_h3k4me3_chipseq'

set :snap_id, `cat SNAPID`.chomp #ec2 eu-west-1 
set :vol_id, `cat VOLUMEID`.chomp #empty until you've created a new volume
set :ebs_size, 20  #Needs to be the size of the snap plus enough space for alignments
set :ebs_zone, 'eu-west-1b'  #is where the ubuntu ami is
set :dev, '/dev/sdf'
set :mount_point, '/mnt/data'

#grab the snapshot ID for the raw data (fastq)
task :get_snap_id, :roles=>:master do
  `curl http://github.com/cassj/ns5dastro_h34kme3_chipseq/raw/master/data/SNAPID > SNAPID `
end 


#make a new EBS volume from this snap 
#cap EBS:create

#and mount your EBS
#cap EBS:attach
#cap EBS:mount_xfs



# There's little point in realigning this data - anything ELAND couldn't align
# has been thrown away already. So




# fetch samtools from svn
desc "get samtools"
task :get_samtools, :roles => :chipseq do
  sudo "apt-get -y install subversion"
  run "svn co https://samtools.svn.sourceforge.net/svnroot/samtools/trunk/samtools"
end
before "get_samtools", "EC2:start"


desc "build samtools"
task :build_samtools, :roles => :chipseq do
  sudo "apt-get -y install zlib1g-dev libncurses5-dev"
  run "cd /home/ubuntu/samtools && make"
end
before "build_samtools", "EC2:start"


desc "install samtools"
task :install_samtools, :roles => :chipseq do
  sudo "cp /home/ubuntu/samtools/samtools /usr/local/bin/samtools"
end
before "install_samtools", "EC2:start"


desc "make bam from sam"
task :to_bam, :roles => :chipseq do
  run "wget -O #{working_dir}/mm9_lengths  'http://github.com/cassj/my_bioinfo_scripts/raw/master/genomes/mm9_lengths'"
  run "samtools view -bt #{working_dir}/mm9_lengths -o #{working_dir}/export.bam #{working_dir}/export.sam"
end
before "to_bam", "EC2:start"

desc "sort bam"
task :sort_bam, :roles => :chipseq do
  run "samtools sort #{working_dir}/export.bam #{working_dir}/export_sorted"
end
before "sort_bam", "EC2:start"

desc "remove duplicates"
task :rmdups, :roles => :chipseq do
  run "samtools rmdup -s #{working_dir}/export_sorted.bam #{working_dir}/export_nodups.bam"
end
before "rmdups", "EC2:start"

desc "index bam files"
task :index, :roles => :chipseq do
  run "samtools index #{working_dir}/export_nodups.bam #{working_dir}/export_nodups.bai"
end
before "index", "EC2:start"



desc "upload to s3"
task "to_s3", :roles => :chipseq do

  servers = find_servers_for_task(current_task)
  puts servers
  it = 0..(servers.length - 1)
  it.each do |i|
    host = servers[i]
    bucket = bucket_names[i]
    object = object_names[i]
    bucket = "bam."+bucket
    run("s3cmd mb s3://#{bucket}", :hosts => host)
    run("s3cmd put #{working_dir}/export.bam s3://#{bucket}", :hosts => host)
    run("s3cmd put #{working_dir}/export_sorted.bam s3://#{bucket}", :hosts => host)
    run("s3cmd put #{working_dir}/export_nodups.bam s3://#{bucket}", :hosts => host)
    run("s3cmd put #{working_dir}/export_nodups.bai s3://#{bucket}", :hosts => host)
  end

end
before "to_s3", "EC2:start"








#if you want to keep the results

#cap EBS:snapshot


#and then shut everything down:

# cap EBS:unmount
# cap EBS:detach
# cap EBS:delete - unless you're planning to use it again.
# cap EC2:stop




