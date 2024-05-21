#!/usr/bin/perl

use strict;
use warnings;

# Folder path as an argument
my $folder_path = $ARGV[0] // die "Usage: $0 folder_path\n";

# Check if the folder path is provided and it is a directory
die "Please provide a directory.\n" if not -d $folder_path;

# Open the directory
opendir my $dh, $folder_path or die "Cannot open directory $folder_path: $!";

# Read all files from the directory excluding '.' and '..'
my @files = grep { -f "$folder_path/$_" } readdir($dh);

# Close the directory handle
closedir $dh;

# Process each file
foreach my $file (@files) {
    # Full path for file
    my $file_path = "$folder_path/$file";

    # Read the file
    open my $in, '<', $file_path or die "Cannot open file $file_path: $!";
    # Create a new file name for the output
    my $out_file_path = "$file_path.filtered";
    open my $out, '>', $out_file_path or die "Cannot open file $out_file_path: $!";

    while (my $line = <$in>) {
        # Check for matching patterns and write to the output file if it matches
        if ($line =~ /Mean:,-?\d+(\.\d+)?([eE][-+]?\d+)?,stdev:,-?\d+(\.\d+)?([eE][-+]?\d+)?,CV:, ?-?\d+(\.\d+)?([eE][-+]?\d+)?,entries:\d+/ ||
            $line =~ /Found cutset, prev_cutset_prefix remains to:,\d+, repetitions:,\d+/ ||
            $line =~ /Iter:,\d+,cutset_reward:,[\d.]+,clean_cutset:,\[(.+)\]/ ||
            $line =~ /phi_\w+:,\d+\.\d+/ ||
            $line =~ /FINISHED CPU \d+\.\d+ MEM \d+ MAXMEM \d+/) {
            print $out $line;
        }
    }

    close $in;
    close $out;

    # Optionally, rename the .filtered file to the original file name to replace it
    rename($out_file_path, $file_path) or die "Cannot rename $out_file_path to $file_path: $!";
}

print "Processing complete. Matched lines have been written to .filtered files\n"

