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
set :ebs_size, 50  #Needs to be the size of the snap plus enough space for alignments
set :ebs_zone, 'eu-west-1b'  #is where the ubuntu ami is
set :dev, '/dev/sdh'
set :mount_point, '/mnt/data2'


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

  #64bit UCSC liftOver
  run "curl http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/liftOver > #{working_dir}/liftOver"
  run "sudo mv #{working_dir}/liftOver /usr/bin/liftOver"
  run 'sudo chmod +x /usr/bin/liftOver'

end 
before "install_liftover", "EC2:start"


desc "get mm8tomm9 chainfile"
task :mm8_chain_mm9, :roles => group_name do
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
    run "cd #{mount_point} && #{working_dir}/scripts/sortedmm8tomm9.R #{f} #{working_dir}/lib/mm8ToMm9.over.chain" 
  }
end
before "liftOver", "EC2:start"


desc "Remove anything mapping to a random chr"
task :remove_random, :roles => group_name do
  files = capture "ls #{mount_point}"
  files = files.split("\n").select{|f| f.match(/unique_hits_mm9\.txt/)}
  files.each{|f|     
    run "cd #{mount_point} && perl -ni.bak -e 'print $_ unless /.*random.*/;' #{f}" 
  }
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
  run "cd #{working_dir}/scripts && curl http://github.com/cassj/ns5dastro_h34kme3_chipseq/raw/master/scripts/sorted2sam.pl > sorted2sam.pl"
  run "sudo mv #{working_dir}/scripts/sorted2sam.pl /usr/local/bin"
  run "sudo chmod +x /usr/local/bin/sorted2sam.pl"
  ip_i = ip.sub('.txt','_mm9.txt')
  input_i = input.sub('.txt','_mm9.txt')
  ip_o = ip_i.sub('.txt','.sam')
  input_o = input_i.sub('.txt','.sam')
  run "sorted2sam.pl #{ip_i} > #{ip_o}"
  run "sorted2sam.pl #{input_i} > #{input_o}"
end
before 'make_sam', 'EC2:start'


desc "add sam headers"
task :sam_head, :roles => group_name do
  ip_i = ip.sub('.txt','_mm9.sam')
  input_i = input.sub('.txt','_mm9.sam')
  ip_o = ip_i.sub('_mm9', '_mm9_head')
  input_o = input_i.sub('_mm9', '_mm9_head')
  upload('scripts/mm9_sam_header', "#{working_dir}/mm9_sam_header")
  run "cd #{mount_point} && cat #{working_dir}/mm9_sam_header #{ip_i} >  #{ip_o}"
  run "cd #{mount_point} && cat #{working_dir}/mm9_sam_header #{input_i} >  #{input_o}"
  run "cd #{mount_point} && mv #{ip_o} #{ip_i}"
  run "cd #{mount_point} && mv #{input_o} #{input_i}"
end
before 'sam_head', 'EC2:start'


desc "make bam from sam"
task :to_bam, :roles => group_name do
  run "cd #{working_dir} && curl http://github.com/cassj/my_bioinfo_scripts/raw/master/genomes/mm9_lengths > mm9_lengths"
  ip_i = ip.sub('.txt','_mm9.sam')
  input_i = input.sub('.txt','_mm9.sam')
  ip_o = ip_i.sub('.sam','.bam')
  input_o = input_i.sub('.sam','.bam')
  run "samtools view -bt #{working_dir}/mm9_lengths -o #{ip_o}  #{ip_i}"
  run "samtools view -bt #{working_dir}/mm9_lengths -o #{input_o} #{input_i}"
end
before "to_bam", "EC2:start"

desc "sort bam"
task :sort_bam, :roles => group_name do
  ip_i = ip.sub('.txt','_mm9.bam')
  input_i = input.sub('.txt','_mm9.bam')
  ip_o = ip_i.sub('.bam','_sorted')
  input_o = input_i.sub('.bam','_sorted')
  run "samtools sort #{ip_i} #{ip_o}"
  run "samtools sort #{input_i} #{input_o}"
end
before "sort_bam", "EC2:start"

desc "remove duplicates"
task :rmdups, :roles => group_name do
  ip_i = ip.sub('.txt','_mm9_sorted.bam')
  input_i = input.sub('.txt','_mm9_sorted.bam')
  ip_o = ip_i.sub('sorted', 'sorted_nodups')
  input_o = input_i.sub('sorted','sorted_nodups')

  run "samtools rmdup -s #{ip_i} #{ip_o}"
  run "samtools rmdup -s #{input_i} #{input_o}"
end
before "rmdups", "EC2:start"


desc "index bam files"
task :index, :roles => group_name do
  ip_i = ip.sub('.txt','_mm9_sorted_nodups.bam')
  input_i = input.sub('.txt','_mm9_sorted_nodups.bam')
  ip_o = ip_i.sub('.bam', '.bai')
  input_o = input_i.sub('.bam','.bai')
  run "samtools index #{ip_i} #{ip_o}"
  run "samtools index #{input_i} #{input_o}"
end
before "index", "EC2:start"


desc "download bam files"
task :get_bam, :roles => group_name do
  `rm -Rf results/alignment/bowtie` #remove previous results
  `mkdir -p results/alignment/bowtie`
  files = capture "ls #{mount_point}"
  files = files.split("\n").select{|f| f.match(/sorted_nodups/)}
  files.each{|f|
    download( "#{mount_point}/#{f}", "results/alignment/bowtie/#{f}")
  }
end
before "get_bam", 'EC2:start'






### Macs ?

macs_url ="http://liulab.dfci.harvard.edu/MACS/src/MACS-1.4.0beta.tar.gz"
macs_version = "MACS-1.4.0beta"

task :install_macs, :roles => group_name do
  sudo "apt-get install -y python"
  run "cd #{working_dir} && wget --http-user macs --http-passwd chipseq #{macs_url}"
  run "cd #{working_dir} && tar -xvzf #{macs_version}.tar.gz"
  run "cd #{working_dir}/#{macs_version} && sudo python setup.py install"
  sudo "ln -s /usr/local/bin/macs* /usr/local/bin/macs"
end
before "install_macs", 'EC2:start'

task :install_peaksplitter, :roles => group_name do
  url ='http://www.ebi.ac.uk/bertone/software/PeakSplitter_Cpp_1.0.tar.gz'
  filename = 'PeakSplitter_Cpp_1.0.tar.gz'
  bin = 'PeakSplitter_Cpp/PeakSplitter_Linux64/PeakSplitter'
  run "cd #{working_dir} && curl #{url} > #{filename}"
  run "cd #{working_dir} && tar -xvzf #{filename}"
  run "sudo cp #{working_dir}/#{bin} /usr/local/bin/PeakSplitter"
end 
before 'install_peaksplitter', 'EC2:start'

#you'll need to have done "install_r" and install_peak_splitter to do this
task :run_macs, :roles => group_name do

  treatment = "#{mount_point}/CMN027_181_unique_hits_mm9_sorted_nodups.bam"
  control = "#{mount_point}/CMN028_028_unique_hits_mm9_sorted_nodups.bam"
  genome = 'mm'
  bws = [300]
  pvalues = [0.00001]

  #unsure what p values and bandwidths are appropriate, try a few?
  bws.each {|bw|
    pvalues.each { |pvalue|

      dir = "#{mount_point}/macs_#{bw}_#{pvalue}"
      run "rm -Rf #{dir}"
      run "mkdir #{dir}"

      macs_cmd =  "macs --treatment #{treatment} --control #{control} --name #{group_name} --format BAM --gsize #{genome} --bw #{bw} --pvalue #{pvalue}"
      run "cd #{dir} && #{macs_cmd}"
      
      dir = "#{mount_point}/macs_#{bw}_#{pvalue}_subpeaks"
      run "rm -Rf #{dir}"
      run "mkdir #{dir}"

      # With SubPeak finding
      # this will take a lot longer as you have to save the wig file 
      macs_cmd =  "macs --treatment #{treatment} --control #{control} --name #{group_name} --format BAM --gsize #{genome} --call-subpeaks  --bw #{bw} --pvalue #{pvalue} --wig"
      run "cd #{dir} && #{macs_cmd}"

    }
  }
  
end
before 'run_macs', 'EC2:start'



#pack up the runs and downloads them to the server (without the wig files)
task :pack_macs, :roles => group_name do
  macs_dirs = capture "ls #{mount_point}"
  macs_dirs = macs_dirs.split("\n").select {|f| f.match(/.*macs.*/)}
  macs_dirs.each{|d|
    run "cd #{mount_point} &&  tar --exclude *_wiggle* -cvzf #{d}.tgz #{d}"
  }
  
end
before 'pack_macs','EC2:start' 

task :get_macs, :roles => group_name do
  macs_files = capture "ls #{mount_point}"
  macs_files = macs_files.split("\n").select {|f| f.match(/.*macs.*\.tgz/)}
  res_dir = 'results/alignment/bowtie/peakfinding/macs'
  `rm -Rf #{res_dir}`
  `mkdir -p #{res_dir}`
  macs_files.each{|f| 
    download("#{mount_point}/#{f}", "#{res_dir}/#{f}") 
    `cd #{res_dir} && tar -xvzf #{f}`
  }

end
before 'get_macs', 'EC2:start'






#if you want to keep the results

#cap EBS:snapshot


#and then shut everything down:

# cap EBS:unmount
# cap EBS:detach
# cap EBS:delete - unless you're planning to use it again.
# cap EC2:stop




