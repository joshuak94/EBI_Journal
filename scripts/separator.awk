BEGIN {
    prev = 1;
    output_pre = "./mappability_separated/mus_musculus";
    output_ext = ".map";
    output_file = output_pre"_1"output_ext;
}

{
    if (prev == $1)
    {
        print $0 >> output_file;
    }
    else
    {
        prev = $1;
        output_file = (output_pre"_"$1output_ext);
        print $0 >> output_file;
    }
}
