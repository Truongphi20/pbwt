import vcfpy


class PBWT:
    def __init__(self, filename):
        self.reader = vcfpy.Reader.from_path(filename)

    def pbwtCursorForwardsA(self):
        # src/pbwtCore.c:460
        pass

    def pbwtReadVcfGT(self):
        # src/pbwtHtslib.c:54
        pass


if __name__ == "__main__":
    filename = "tests/data/OMNI.vcf"
    
    pbwt = PBWT(filename)