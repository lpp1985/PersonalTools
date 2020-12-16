# 导入 -> 系统 package
use strict;
use warnings;
use File::Spec;
use Getopt::Long;

# 定义 -> 常量
use constant SCRIPTDIR => (File::Spec->splitpath(File::Spec->rel2abs($0)))[1];
use constant PWD => $ENV{"PWD"};

# 定义 -> 核心变量


# 检测 -> 脚本输入
my ($gdc_sample_sheet, $download_dir, $output_file, $gene_info, $if_help);
GetOptions(
    "gdc_sample_sheet|s=s"     => \$gdc_sample_sheet,
    "download_dir|d=s"         => \$download_dir,
    "output_file|o=s"          => \$output_file,
    "gene_info|g=s"            => \$gene_info,
    "help|h"                   => \$if_help,
);
die help() if(defined $if_help or (not defined $gdc_sample_sheet or not defined $download_dir or not defined $output_file));
###################################################################### 主程序

# 读入meta信息，以样本名作为第一键值
my %hashMetadata = read_matrix($gdc_sample_sheet, 0);
check_duplicate_sample_id(\%hashMetadata);

# 读入gene_info 信息，用于注释
my %hashGeneInfo;
my @anno_heads;
if(defined $gene_info and is_file_ok($gene_info) == 1)
{   
    %hashGeneInfo = read_matrix($gene_info, 0);
    @anno_heads   = @{$hashGeneInfo{'Head'}};
}
 


# 读入表达量
my @ids = sort keys %{$hashMetadata{'Data'}};
my @samples;
my %hashExpression;
print "读取表达量数据\n";
foreach my $id( @ids )
{   
    process_bar_array($id, \@ids);
    my $file_id   = $hashMetadata{'Data'}{$id}{'File ID'};
    my $file_name = $hashMetadata{'Data'}{$id}{'File Name'};
    my $sample_id = $hashMetadata{'Data'}{$id}{'Sample ID'};
    push @samples, $sample_id;

    my $file = "$download_dir/$file_id/$file_name";
    die "[Error] 文件丢失：$file" if(is_file_ok($file) == 0);

    open FILE, "gzip -cd $file|";
    while(<FILE>)
    {
        $_=~s/[\r\n]//g;
        next if($_!~/\w/);
        next if($_=~/^__/);  # 去掉特殊记录 __no_feature __ambiguous __too_low_aQual __not_aligned __alignment_not_unique
        my ($gene, $value) = split /\t/, $_;
        $hashExpression{$gene}{$sample_id} = $value;
    }
    close FILE;
}

# 输出
print "output $output_file\n";
open OUTPUT, ">$output_file";
print OUTPUT "Gene\t" . (join "\t", (@samples, @anno_heads)) . "\n";
foreach my $gene( sort keys %hashExpression)
{
    my @exps = map{ $hashExpression{$gene}{$_} } @samples;
    my @annos = ();
       @annos = map{ $hashGeneInfo{'Data'}{$gene}{$_} } @anno_heads if(length(@anno_heads) > 0 );
    
    print OUTPUT "$gene\t" . (join "\t", (@exps, @annos)) . "\n";
}
close OUTPUT;


###################################################################### 子函数

# 检验文件是否为空
sub is_file_ok{
    my @files = @_;
    my $isOK = 1;
    foreach my $file(@files)
    {
        $isOK = 0 if (not -e $file or not -f $file or -s $file == 0);
    }    
    return $isOK;
}

# 读取带有表头的矩阵
sub read_matrix{
    my $file = shift @_;
    my $key_col = shift @_;

    print "Read $file\n";
    my %hashMatrix;
    open FILE, $file;
    my $line1 = <FILE>;
       $line1 =~ s/[\r\n]//g;
    my @heads = split /\t/, $line1;
    $hashMatrix{'Head'} = \@heads;

    while(<FILE>)
    {
        $_ =~ s/[\r\n]//g;
        my @datas = split /\t/, $_;
        my $key_value = $datas[$key_col];
        foreach my $col(0..$#heads)
        {
            my $value = exists $datas[$col] ? $datas[$col] : "";
            $hashMatrix{'Data'}{$key_value}{$heads[$col]} = $value;
        }
    }
    close FILE;
    return %hashMatrix;
}

# 进度条， 根据输入向量计算
sub process_bar_array{
    my $process_name   = shift @_; # 需要分析的对象
    my $process_arrays = shift @_; # 总列表
    my $process_all_count = @$process_arrays;# 需要分析的总列表数量
    my ($pos) = grep $$process_arrays[$_] eq $process_name, 0..$process_all_count-1;
    my $process_count = $pos+ 1; # 当前对象的位置
    my $process_perc = sprintf "%0.2f", 100 * $process_count/$process_all_count; # 进度
    my $windows = 100; # 窗口宽度
    my $finished = $windows * int($process_perc) / 100; # 完成
    my $unfinished = $windows - $finished; # 未完成
    print "\r[". '>' x $finished . ' ' x $unfinished . "]  [$process_count/$process_all_count]  $process_perc% ";
    print "\n" if($process_count == $process_all_count);   
}
sub check_duplicate_sample_id{
    my $hashMetadata = shift @_;
    my %hashSample;
    foreach my $id(keys %{$hashMetadata->{'Data'}})
    {
        my $sample_id = $hashMetadata->{'Data'}{$id}{'Sample ID'};
        $hashSample{$sample_id}++;
    }
    my @dups = grep{ $hashSample{$_} > 1  } sort keys %hashSample;
    die "[Error] gdc_sample_sheet 文件中有重复的 Sample ID, 请确保唯一： @dups\n" if(scalar(@dups) > 0 );
}

sub help{
    my $info = "
Program: 合并tcga下载的fpkm文件
Version: 2020-08-05
Contact: 129 甘斌

Usage:   perl ".(File::Spec->splitpath(File::Spec->rel2abs($0)))[2]." [options]

Options: [必填]

         --gdc_sample_sheet/-s      在TCGA的‘Cart’中，下载的文件sample sheet。重点需要三个表头信息：File ID，File Name，Sample ID。 文件的第一列务必保证不要有重复，因为代码是用第一列作为哈希的键值。
                                    注意检查Sample ID 是否有重复。
         --download_dir/-d          gdc-client下载文件的保存路径
         --output_file/-o           输出文件

Options: [选填]
         --gene_info/-g     基因注释文件，可通过TCGA下载https://api.gdc.cancer.gov/data/b011ee3e-14d8-4a97-aed4-e0b10f6bbe82
                            官方注释文件：gencode.gene.info.v22.tsv
         --help/-h          查看帮助文档
    \n";
    return $info;
}

