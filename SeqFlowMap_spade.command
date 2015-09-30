#!/bin/sh

dir=${0%/*}
if [ -d "$dir" ]; then
    cd "$dir"
fi

sed -e 's/panel_files.*/panel_files=dir(pattern=glob2rx("*.fcs")),/' -e 's/LAYOUT_TABLE.*//' -e 's/MST_GRAPH.*//' -e 's/SPADE.plot.trees.*//' <runSPADE.R >tmp.R

for file in *.fcs
do
    mkdir "$file"_dir
    mv "$file" "$file"_dir/"$file"
    cp tmp.R "$file"_dir/runSPADE_mod.R
    cd "$file"_dir/
    R -f runSPADE_mod.R
    cd ..
done

rm tmp.R