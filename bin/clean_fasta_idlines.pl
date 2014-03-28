#!/usr/bin/perl -w
use strict;

while(<>){
s/^(>\s*\S+).*/$1/; # just keep > and id.
print;
}
