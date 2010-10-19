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
set :ebs_size, 30  #Needs to be the size of the snap plus enough space for alignments
set :ebs_zone, 'eu-west-1b'  #is where the ubuntu ami is
set :dev, '/dev/sdf'
set :mount_point, '/mnt/data'


set :ip, "#{mount_point}/CMN027_181_unique_hits.txt"
set :input, "#{mount_point}/CMN028_028_unique_hits.txt"


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
# has been thrown away already. So...

desc "remove the .fa extension from chr names in unique_hits files"
task :strip_fa, :roles => group_name do
  run "perl -pi -e 's/(chr.+)\.fa/$1/g' #{ip}\n"
  run "perl -pi -e 's/(chr.+)\.fa/$1/g' #{input}\n"

end
before "strip_fa", "EC2:start"



desc "install R on all running instances in group group_name"
task :install_r, :roles  => group_name do
  sudo "mkdir -p #{working_dir}/scripts"
  user = variables[:ssh_options][:user]
  sudo "chown #{user} #{working_dir}/scripts"
  sudo 'apt-get -y install r-base'
  sudo 'apt-get -y install build-essential libxml2 libxml2-dev libcurl3 libcurl4-openssl-dev'
  run "curl http://github.com/cassj/ns5dastro_h34kme3_chipseq/raw/master/scripts/R_setup.R > #{working_dir}/scripts/R_setup.R"
  sudo "Rscript #{working_dir}/scripts/R_setup.R"
end
before "install_r", "EC2:start"
  


desc "install liftOver"
task :install_liftover, :roles => group_name do

  #liftover
  sudo 'wget -O /usr/bin/liftOver http://hgdownload.cse.ucsc.edu/admin/exe/linux.i386/liftOver'
  sudo 'chmod +x /usr/bin/liftOver'

end 
before "install_liftover", "EC2:start"

desc "get mm8tomm9 chainfile"
task :mm8_chain_mm9, :roles => group_name do
  #chainfile - note that the script sortedmm8tomm9 expects the chain files to be in ../lib
  run "mkdir -p #{working_dir}/lib"
  run "curl http://hgdownload.cse.ucsc.edu/goldenPath/mm8/liftOver/mm8ToMm9.over.chain.gz > #{working_dir}/lib/mm8ToMm9.over.chain.gz"
  run "gunzip -c #{working_dir}/lib/mm8ToMm9.over.chain.gz > #{working_dir}/lib/mm8ToMm9.over.chain"
end 
before "mm8_chain_mm9", "EC2:start"



desc "liftOver ELAND mm8 positions to mm9"
task :liftOver, :roles => group_name do

  run "curl  http://github.com/cassj/ns5dastro_h34kme3_chipseq/raw/master/scripts/sortedmm8tomm9.R > #{working_dir}/scripts/sortedmm8tomm9.R"
  run "chmod +x #{working_dir}/scripts/sortedmm8tomm9.R"
 
  files = capture "ls #{mount_point}"
  files = files.split("\n").select{|f| f.match(/unique_hits\.txt/)}
  files.each{|f|     
    run "cd #{mount_point} && #{working_dir}/scripts/sortedmm8tomm9.R #{f}" 
  }
end
before "liftOver", "EC2:start"

#I am here

desc "Remove anything mapping to a random chr"
task :remove_random, :roles => group_name do
  run "cd #{working_dir} && perl -ni.bak -e 'print $_ unless /.*random.*/;' *_mm9.txt"
end
before "remove_random", "EC2:start"




# fetch samtools from svn
desc "get samtools"
task :get_samtools, :roles => group_name do
  sudo "apt-get -y install subversion"
  run "svn co https://samtools.svn.sourceforge.net/svnroot/samtools/trunk/samtools"
end
before "get_samtools", "EC2:start"


desc "build samtools"
task :build_samtools, :roles => group_name do
  sudo "apt-get -y install zlib1g-dev libncurses5-dev"
  run "cd /home/ubuntu/samtools && make"
end
before "build_samtools", "EC2:start"


desc "install samtools"
task :install_samtools, :roles => group_name do
  sudo "cp /home/ubuntu/samtools/samtools /usr/local/bin/samtools"
end
before "install_samtools", "EC2:start"



desc "make sam files from unique.txt"
task :make_sam, :roles => group_name do
  run "cd #{working_dir} && curl http://github.com/cassj/ns5dastro_h34kme3_chipseq/raw/master/scripts/sorted2sam.pl > sorted2sam.pl"
  run "sudo mv #{working_dir}/sorted2sam.pl /usr/local/bin"
  run "sudo chmod +x /usr/local/bin/sorted2sam.pl"
  ip_o = ip.sub('.txt','.sam')
  input_o = input.sub('.txt','.sam')
  run "sorted2sam.pl #{ip} > #{ip_o}"
  run "sorted2sam.pl #{input} > #{input_o}"
end
before 'make_sam', 'EC2:start'


desc "make bam from sam"
task :to_bam, :roles => group_name do
  run "cd #{working_dir} && curl http://github.com/cassj/my_bioinfo_scripts/raw/master/genomes/mm9_lengths > mm9_lengths"
  ip_i = ip.sub('.txt','.sam')
  input_i = input.sub('.txt','.sam')
  ip_o = ip.sub('.txt','.bam')
  input_o = input.sub('.txt','.bam')
  run "samtools view -bt #{working_dir}/mm9_lengths -o #{ip_o}  #{ip_i}"
  run "samtools view -bt #{working_dir}/mm9_lengths -o #{input_o} #{input_i}"
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




