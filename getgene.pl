#!/usr/bin/perl -w
use strict;
use List::MoreUtils ':all';
use List::Util qw(min max);

# 初始化数据结构
my (%chr,  %pos, @genes4spall, @genes4spshare);
my ($i, $j, $gene, $file, $out, $line, $path);
my @col = ();
my @tmp = ();
my %overlap_genes=();
# 读取染色体长度文件
for $i ("anak", "blochii", "carolinus", "ovatus") {
    $file = $i . ".chr.length";
    open IN, $file or die "无法打开 $file: $!";
    while ($line = <IN>) {
        chomp $line;
        @col = split /\t/, $line;
        $chr{$i}{$col[0]} = $col[1];
    }
    close IN;
    print STDERR "完成读取 $file\n";
}

# 处理每个物种的完整表格
for $i ("anak", "blochii", "carolinus", "ovatus") {
    $file = $i . ".full_table.tsv";
    open IN, $file or die "无法打开 $file: $!";
    
    # 存储所有基因区间
    my @all_intervals;
    my %gene_to_interval;
    
    while ($line = <IN>) {
        chomp $line;
        next if $line =~ /^#/;
        if ($line =~ /Complete/) {
            @col = split /\t/, $line;
            if ($col[3] >= 1000 && $col[3] <= ($chr{$i}{$col[2]} - 1000) &&
                $col[4] >= 1000 && $col[4] <= ($chr{$i}{$col[2]} - 1000)) {
                
                my ($start, $end) = ($col[3] < $col[4]) ? ($col[3], $col[4]) : ($col[4], $col[3]);
                
                # 存储区间信息
                my $interval = {
                    chr   => $col[2],
                    start => $start,
                    end   => $end,
                    gene  => $col[0]
                };
                
                push @all_intervals, $interval;
                $gene_to_interval{$col[0]} = $interval;
                
                push @genes4spall, $col[0];
            }
        }
    }
    close IN;
    
    # 检测重叠区间
    
    # 按染色体分组
    my %intervals_by_chr;
    foreach my $interval (@all_intervals) {
        push @{$intervals_by_chr{$interval->{chr}}}, $interval;
    }
    
    # 对每个染色体上的区间排序
    foreach my $chr (keys %intervals_by_chr) {
        # 按起始位置排序
        my @sorted_intervals = sort { $a->{start} <=> $b->{start} } @{$intervals_by_chr{$chr}};
        
        # 检测重叠
        for (my $k = 0; $k < @sorted_intervals; $k++) {
            my $current = $sorted_intervals[$k];
            
            # 检查与后续区间的重叠
            for (my $m = $k + 1; $m < @sorted_intervals; $m++) {
                my $next = $sorted_intervals[$m];
                
                # 如果下一个区间在当前区间之后开始，则停止检查
                last if $next->{start} > $current->{end};
                
                # 检测到重叠
                if ($next->{start} < $current->{end}) {
                    $overlap_genes{$current->{gene}} = 1;
                    $overlap_genes{$next->{gene}} = 1;
                }
            }
        }
    }
    
    
    foreach my $interval (@all_intervals) {
        if(exists $overlap_genes{$interval->{gene}}){
		next;
	}
	else{
        	$pos{$i}{$interval->{gene}}=$interval->{chr}."\t".$interval->{start}."\t".$interval->{end};
        }
    }
    
    my $total_genes = scalar @all_intervals;
    my $overlap_count = scalar keys %overlap_genes;
    print STDERR "物种 $i: 总基因数 = $total_genes, 重叠基因数 = $overlap_count\n";
}

# 识别共享基因
@genes4spall = uniq @genes4spall;
foreach $i (@genes4spall) {
    if(defined($pos{"anak"}{$i}) && 
        defined($pos{"blochii"}{$i}) && 
        defined($pos{"carolinus"}{$i}) && 
        defined($pos{"ovatus"}{$i})){
        	push @genes4spshare, $i;
	}
}
print STDERR "共享基因数量: ", scalar @genes4spshare, "\n";

# 输出每个物种的共享基因位置
foreach $i ("anak", "blochii", "carolinus", "ovatus") {
    $out = $i . ".txt";
    open OUT, ">$out" or die "无法创建 $out: $!";
    foreach $j (@genes4spshare) {
        print OUT $pos{$i}{$j}, "\n" if defined $pos{$i}{$j};
    }
    close OUT;

    # 输出4D位点
    $out = $i . ".4Dsites.txt";
    open OUT, ">$out" or die "无法创建 $out: $!";
    foreach $j (@genes4spshare) {
        $file = $i . ".gff/" . $j . ".4Dsite.gz";
        open IN, "gunzip -cd $file |" or warn "无法打开 $file: $!";
        while ($line = <IN>) {
            chomp $line;
            @col = split /\t/, $line;
            print OUT $col[0], "\t", $col[1], "\n";
        }
        close IN;
    }
    close OUT;
}
