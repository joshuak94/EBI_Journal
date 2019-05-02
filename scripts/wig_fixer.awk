BEGIN {
    out_file = "filtered_mus_musculus_all_chrom.map.wig";
}
{
    chr = $1;
    if ($1 == "fixedStep")
    {
        chr = substr($2, 7);
        if (chr == "1" || chr == "2" || chr == "3" || chr == "4" || chr == "5" || chr == "6" || chr == "7" || chr == "8" || chr == "9" || chr == "10" || chr == "11" || chr == "12" || chr == "13" || chr == "14" || chr == "15" || chr == "16" || chr == "17" || chr == "18" || chr == "19" || chr == "X" || chr == "Y")
        {
            prevStart = substr($3, 7);
            prevStep = substr($4, 6);
            prevChr = substr($2, 7);
        }
    }
    else if ($1 == "1.000000")
    {
        print (prevChr"\t"(prevStart)"\t"(prevStart + prevStep)"\t1") >> out_file;
        prevStart = prevStart + prevStep;
    }
    else if (chr == "1" || chr == "2" || chr == "3" || chr == "4" || chr == "5" || chr == "6" || chr == "7" || chr == "8" || chr == "9" || chr == "10" || chr == "11" || chr == "12" || chr == "13" || chr == "14" || chr == "15" || chr == "16" || chr == "17" || chr == "18" || chr == "19" || chr == "X" || chr == "Y")
    {
        print $0 >> out_file;
    }
}
