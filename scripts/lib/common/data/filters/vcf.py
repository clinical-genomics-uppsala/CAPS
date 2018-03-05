import vcf.filters

class AlleleRationFilterDP(vcf.filters.Base):
    "Filter sites by allele_ratio using DP field"

    name = 'ad'

    @classmethod
    def customize_parser(self, parser):
        parser.add_argument('--min-allele-ratio', type=float, default=0.01,
                help='Filter sites below this allele_ratio')

    def __init__(self, args):
        self.min_allele_ratio = args.min_allele_ratio

    def __call__(self, record):
        if record. < self.min_allele_ratio:
            return record.QUAL

class AlleleRationFilterAD(vcf.filters.Base):
    "Filter sites by allele_ratio"

    name = 'ad'

    @classmethod
    def customize_parser(self, parser):
        parser.add_argument('--min-allele-ratio', type=float, default=0.01,
                help='Filter sites below this allele_ratio using AD field')

    def __init__(self, args):
        self.min_allele_ratio = args.min_allele_ratio

    def __call__(self, record):
        record.info['AD'])
        if record. < self.min_allele_ratio:
            return record.QUAL
