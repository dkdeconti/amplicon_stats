def create_report(snps, indels, snp_vcfs, indel_vcfs, bams, plots, dir_map,
                  dry=False):
    '''
    Injects data into html template with jinja2.
    '''
    this_dir = os.path.dirname(os.path.realpath(__file__))
    lib_dir = os.path.join(this_dir, 'lib')
    report_dir = dir_map["reportdir"]
    lib_destination = os.path.join(report_dir, 'lib')
    report = '/'.join([report_dir, "html_report.html"])
    env = jinja2.Environment(loader=jinja2.FileSystemLoader(this_dir))
    template = env.get_template("template.html")
    for vcf, samplename in snp_vcfs.items():
        print vcf, samplename
    samples = {samplename : {"header" : snps[samplename][0],
                             "variants": snps[samplename][1:] + \
                                         indels[samplename][1:],
                             "plots": [p.split('/')[-1]
                                       for p in plots[samplename]],
                             "samplename": samplename.split('/')[-1],
                             "bamfile": os.path.relpath(bams[samplename],
                                                        dir_map["reportdir"]),
                             "snpvcf": os.path.relpath(snp_vcfs[samplename],
                                                       dir_map["reportdir"]),
                             "indelvcf": os.path.relpath(indel_vcfs[samplename],
                                                         dir_map["reportdir"])}
               for samplename in snps.keys()}
    context = {"samples": samples}
    if not dry:
        with open(report, 'w') as outfile:
            outfile.write(template.render(context))
        shutil.copytree(lib_dir, lib_destination)