#!/usr/bin/perl
use strict;
use warnings;

# 声明所有全局变量
my %target_positions;
my %total_count;
my %het_count;
my @samples;

# 参数检查
die "Usage: $0 <vcf.gz> <positions.txt>\n" unless @ARGV == 2;
my ($vcf_gz, $positions_file) = @ARGV;

# 读取目标位点
open my $pos_fh, '<', $positions_file or die "Cannot open $positions_file: $!";
while (<$pos_fh>) {
    chomp;
    my ($chrom, $pos) = split /\t/;
    $target_positions{"$chrom:$pos"} = 1 if defined $chrom && defined $pos;
}
close $pos_fh;

# 错误日志
open my $err_fh, '>', 'vcf_parse_errors.log' or die "Cannot open error log: $!";
print $err_fh "Line\tError\n";

# 处理 VCF 文件
open IN,"gunzip -cd  $vcf_gz |" or die "Can't open $vcf_gz\n";

my $line_count = 0;

while (my $line=<IN>) {
    chomp($line);
    $line_count++;
    
    # 跳过注释行
    if ($line=~/^##/) {
	    next;
    }
    elsif($line=~/^#CHROM/){   # 解析样本名
            my @fields = split(/\t/,$line);
            @samples = @fields[ 9.. $#fields] if @fields > 9;
            # 初始化计数器
            $total_count{$_} = 0 for @samples;
            $het_count{$_} = 0 for @samples;
    }
    else{ 
    	# 安全解析数据行
   	my @fields = split(/\t/,$line);
    
   	# 跳过字段不足的行
    	if (@fields < 9) {
      		print $err_fh "$line_count\tInsufficient fields (".scalar(@fields).")\n";
        	next;
    	}
        next if ($line=~/INDEL/); 
    	my ($chrom, $pos, $id, $ref, $alt, $qual, $filter, $info, $format) = @fields[0..8];
    	my @genotypes = @fields[9..$#fields];
    
    	# 检查变量定义
   	 unless (defined $ref && defined $alt) {
        	print $err_fh "$line_count\tUndefined ref/alt\n";
        	next;
   	 }
    
    	# 只处理 SNP
    	next unless length($ref) == 1 && length($alt) == 1;
    
    	# 只处理目标位点
    	next unless exists $target_positions{"$chrom:$pos"};
    
    	# 解析格式字段
    	my @format_fields = split /:/, $format;
    	my ($gt_index, $dp_index) = (-1, -1);
    	for my $i (0 .. $#format_fields) {
        	$gt_index = $i if $format_fields[$i] eq 'GT';
        	$dp_index = $i if $format_fields[$i] eq 'DP';
    	}
    
    	# 处理每个样本
    	for my $i (0 .. $#samples) {
        	next unless $i <= $#genotypes;  # 防止数组越界
        
        	my $sample = $samples[$i];
        	my @sample_fields = split /:/, $genotypes[$i];
        
        	# 跳过缺失基因型
        	next if $gt_index < 0 || !defined $sample_fields[$gt_index] || $sample_fields[$gt_index] =~ /\./;
        
        	# 检查覆盖度
        	if ($dp_index >= 0 && defined $sample_fields[$dp_index] && $sample_fields[$dp_index] >= 5) {
            		$total_count{$sample}++;
            
            		# 检查杂合性
            		if ($sample_fields[$gt_index] =~ /^(0\/1|1\/0|0\|1|1\|0)$/) {
                		$het_count{$sample}++;
            		}
        	}
    	}
    }
}

close IN;
print $err_fh "The total line number: ",$line_count,"\n";
close $err_fh;
# 输出结果
print "Sample\tTotal_Sites\tHeterozygous_Sites\tHeterozygosity\n";
for my $sample (@samples) {
    my $total = $total_count{$sample} || 0;
    my $het = $het_count{$sample} || 0;
    my $het_ratio = $total > 0 ? sprintf("%.5f", $het/$total) : "NA";
    
    print join("\t", $sample, $total, $het, $het_ratio), "\n";
}
