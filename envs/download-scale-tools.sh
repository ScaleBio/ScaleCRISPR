#!/bin/sh
srcDir=$(dirname -- "$0")
dir=${1:-"$(readlink -f "${srcDir}/../bin")"}
mkdir -p $dir
echo "Downloading to $dir"
curl http://scale.pub.s3.amazonaws.com/scaleTools/cicd/bcParser/master/250403-g84327ca/bc_parser -o $dir/bc_parser && chmod 755 $dir/bc_parser
curl http://scale.pub.s3.amazonaws.com/scaleTools/cicd/seqTagFinder/remove_take/250806-g589a44d/seq_tag_finder -o $dir/seq_tag_finder && chmod 755 $dir/seq_tag_finder
