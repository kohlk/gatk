package org.broadinstitute.hellbender.tools.spark.sv.discovery.inference;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.TextCigarCodec;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.SvType;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.TestUtilsForAssemblyBasedSVDiscovery;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AlignmentInterval;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.ContigAlignmentsModifier;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.testng.Assert;
import org.testng.annotations.Test;
import scala.Tuple2;

import java.util.Collections;
import java.util.List;

import static org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AssemblyContigWithFineTunedAlignments.NO_GOOD_MAPPING_TO_NON_CANONICAL_CHROMOSOME;

public class AssemblyBasedSVDiscoveryTestDataProviderForBreakEndVariants {

    // TODO: 5/23/18 use to fill up AssemblyBasedSVDiscoveryTestDataProviderForBreakEndVariants
    @Test(groups = "sv")
    public void testRefOrderSwitch() {
        AlignmentInterval region1 = new AlignmentInterval(
                // assigned from chr18 to chr21 to use the dict
                new SimpleInterval("chr21", 39477098, 39477363),
                1 ,268,
                TextCigarCodec.decode("236M2I30M108S"), true, 32, 25, 133, ContigAlignmentsModifier.AlnModType.NONE);
        AlignmentInterval region2 = new AlignmentInterval(
                new SimpleInterval("chr21", 39192594, 39192692),
                252 ,350,
                TextCigarCodec.decode("251S99M26S"), true, 32, 1, 94, ContigAlignmentsModifier.AlnModType.NONE);
        SimpleChimera simpleChimera = new SimpleChimera(region1, region2, Collections.emptyList(), "testContig", NO_GOOD_MAPPING_TO_NON_CANONICAL_CHROMOSOME, TestUtilsForAssemblyBasedSVDiscovery.b38_seqDict_chr20_chr21);
        NovelAdjacencyAndAltHaplotype breakpoints = new NovelAdjacencyAndAltHaplotype(simpleChimera,
                "TTCCTTAAAATGCAGGTGAATACAAGAATTAGGTTTCAGGTTTTATATATATATTCTGATATATATATATAATATAACCTGAGATATATATATAAATATATATATTAATATATATTAATATATATAAATATATATATATTAATATATATTTATATATAAATATATATATATTAATATATATAAATATATATAAATATATATATATTAATATATATTAATATATAAATATATATATATTAATATATATTAATATATATAAATATATATATTAATATATATAAATATATATATAAATATATATAAATATATAAATATATATATAAATATATATAAATATATATAAATATATATACACACATACATACACATATACATT".getBytes(),
                TestUtilsForAssemblyBasedSVDiscovery.b38_seqDict_chr20_chr21);
        Assert.assertEquals(breakpoints.getLeftJustifiedLeftRefLoc(), new SimpleInterval("chr21", 39192594, 39192594));
        Assert.assertEquals(breakpoints.getLeftJustifiedRightRefLoc(), new SimpleInterval("chr21", 39477346, 39477346));
        Assert.assertEquals(breakpoints.getComplication().getHomologyForwardStrandRep(), "ATATATAAATATATATA");
        Assert.assertTrue(breakpoints.getComplication().getInsertedSequenceForwardStrandRep().isEmpty());
    }

    public static final class TestDataBreakEndVariants extends AssemblyBasedSVDiscoveryTestDataProvider.AssemblyBasedSVDiscoveryTestDataForSimpleChimera {

        private TestDataBreakEndVariants(final AlignmentInterval firstAlignment, final AlignmentInterval secondAlignment,
                                         final String evidenceAssemblyContigName, final byte[] evidenceContigSeq,
                                         final boolean expectedFirstContigRegionHasLaterReferenceMapping,
                                         final SimpleChimera manuallyCuratedSimpleChimera,
                                         final NovelAdjacencyAndAltHaplotype manuallyCuratedBiPathBubble,
                                         final List<SvType> manuallyCuratedSVTypes,
                                         final List<VariantContext> manuallyCuratedVariants,
                                         final Class<? extends BreakpointsInference> inferencer) {
            super(firstAlignment, secondAlignment, evidenceAssemblyContigName, evidenceContigSeq, expectedFirstContigRegionHasLaterReferenceMapping, manuallyCuratedSimpleChimera, manuallyCuratedBiPathBubble, manuallyCuratedSVTypes, manuallyCuratedVariants, inferencer);
        }

        @Override
        public SAMSequenceDictionary getAppropriateDictionary() {
            return TestUtilsForAssemblyBasedSVDiscovery.bareBoneHg38SAMSeqDict;
        }

        @Override
        public Class<? extends BreakpointsInference> getAppropriateBreakpointInferencer() {
            return inferencer;
        }
    }

    public static final boolean testDataInitialized;

//    public static final TestDataBreakEndVariants forIntraChromosomeRefOrderSwap_plus;
//    public static final TestDataBreakEndVariants forIntraChromosomeRefOrderSwap_minus;
//
//    public static final TestDataBreakEndVariants forIntraChromosomeStrandSwitch_plus;
//    public static final TestDataBreakEndVariants forIntraChromosomeStrandSwitch_minus;
//
//    public static final TestDataBreakEndVariants forInterChromosomeSimple_plus;
//    public static final TestDataBreakEndVariants forInterChromosomeSimple_minus;
//
//    public static final TestDataBreakEndVariants forInterChromosomeStrandSwitch_plus;
//    public static final TestDataBreakEndVariants forInterChromosomeStrandSwitch_minus;

    static {
        testDataInitialized = true;
    }

    public static List<TestDataBreakEndVariants> getAllTestData() {
        return Collections.emptyList();
    }

    public static List<Tuple2<TestDataBreakEndVariants, TestDataBreakEndVariants>> getAllTestDataPaired() {
        return Collections.emptyList();
    }

    //        // same-chr translocation suspect, forward and reverse representation
//        AlignmentInterval intervalOne = new AlignmentInterval(new SimpleInterval("chr20", 61015129, 61015272), 1, 144, TextCigarCodec.decode("144M148H"), true, 60, 1, 139, ContigAlignmentsModifier.AlnModType.NONE);
//        AlignmentInterval intervalTwo = new AlignmentInterval(new SimpleInterval("chr20", 60992732, 60992880), 144, 292, TextCigarCodec.decode("143S149M"), true, 60, 0, 149, ContigAlignmentsModifier.AlnModType.NONE);
//        result.add(new TestData(intervalOne, intervalTwo, bareBoneHg38SAMSeqDict, testData.distances, true));
//
//        intervalOne = new AlignmentInterval(new SimpleInterval("chr20", 28861368, 28861775), 1, 409, TextCigarCodec.decode("387M1I21M623H"), false, 60, 22, 286, ContigAlignmentsModifier.AlnModType.NONE);
//        intervalTwo = new AlignmentInterval(new SimpleInterval("chr20", 28896473, 28897229), 276, 1032, TextCigarCodec.decode("275S757M"), false, 60, 1, 752, ContigAlignmentsModifier.AlnModType.NONE);
//        result.add(new TestData(intervalOne, intervalTwo, bareBoneHg38SAMSeqDict, testData.distances, true));
//
//        // diff-chr translocation suspect without SS
//        intervalOne = new AlignmentInterval(new SimpleInterval("chr21", 24923683, 24923715), 1, 33, TextCigarCodec.decode("33M130H"), true, 60, 0, 33, ContigAlignmentsModifier.AlnModType.NONE);
//        intervalTwo = new AlignmentInterval(new SimpleInterval("chr20", 11590055, 11590197), 21, 163, TextCigarCodec.decode("20S143M"), true, 60, 3, 128, ContigAlignmentsModifier.AlnModType.NONE);
//        result.add(new TestData(intervalOne, intervalTwo, bareBoneHg38SAMSeqDict, testData.distances, true));
//
//        // diff-chr translocation suspect with SS
//        intervalOne = new AlignmentInterval(new SimpleInterval("chr21", 5374092, 5374747), 1, 656, TextCigarCodec.decode("656M322S"), true, 60, 14, 586, ContigAlignmentsModifier.AlnModType.NONE);
//        intervalTwo = new AlignmentInterval(new SimpleInterval("chr20", 28764673, 28765145), 506, 978, TextCigarCodec.decode("473M505H"), false, 60, 16, 393, ContigAlignmentsModifier.AlnModType.NONE);
//        result.add(new TestData(intervalOne, intervalTwo, bareBoneHg38SAMSeqDict, testData.distances, true));
//
//        // same-chr reference order switch, but overlaps (hence incomplete picture)
//        intervalOne = new AlignmentInterval(new SimpleInterval("20", 283, 651), 383, 751, TextCigarCodec.decode("382H369M274H"), true, 60, 23, 254, ContigAlignmentsModifier.AlnModType.NONE);
//        intervalTwo = new AlignmentInterval(new SimpleInterval("20", 1, 413), 613, 1025, TextCigarCodec.decode("612H413M"), true, 60, 0, 413, ContigAlignmentsModifier.AlnModType.NONE);
//        result.add(new TestData(intervalOne, intervalTwo, b37_seqDict, testData.distances, true));
//        // same-chr translocation suspect, forward and reverse representation
//        data.add(new Object[]{tuple3s.get(i), StrandSwitch.NO_SWITCH, true, new Tuple2<>(tuple3s.get(i)._2().referenceSpan, tuple3s.get(i)._1().referenceSpan)}); ++i;
//
//        data.add(new Object[]{tuple3s.get(i), StrandSwitch.NO_SWITCH, false, new Tuple2<>(tuple3s.get(i)._1().referenceSpan, tuple3s.get(i)._2().referenceSpan)}); ++i;
//
//        // diff-chr translocation suspect without SS
//        data.add(new Object[]{tuple3s.get(i), StrandSwitch.NO_SWITCH, true, new Tuple2<>(tuple3s.get(i)._2().referenceSpan, tuple3s.get(i)._1().referenceSpan)}); ++i;
//
//        // diff-chr translocation suspect with SS
//        data.add(new Object[]{tuple3s.get(i), StrandSwitch.FORWARD_TO_REVERSE, false, new Tuple2<>(tuple3s.get(i)._2().referenceSpan, tuple3s.get(i)._1().referenceSpan)}); ++i;
//
//        // same-chr reference order switch, but overlaps (hence incomplete picture)
//        data.add(new Object[]{tuple3s.get(i), StrandSwitch.NO_SWITCH, true, new Tuple2<>(tuple3s.get(i)._2().referenceSpan, tuple3s.get(i)._1().referenceSpan)}); ++i;
//    private static Tuple2<TestDataBreakEndVariants, TestDataBreakEndVariants> forIntraChromosomeRefOrderSwap_plus() {
//
//    }
//    private static Tuple2<TestDataBreakEndVariants, TestDataBreakEndVariants> forIntraChromosomeRefOrderSwap_minus() {
//
//    }
//
//
//    private static Tuple2<TestDataBreakEndVariants, TestDataBreakEndVariants> forIntraChromosomeStrandSwitch_plus() {
//
//    }
//    private static Tuple2<TestDataBreakEndVariants, TestDataBreakEndVariants> forIntraChromosomeStrandSwitch_minus() {
//
//    }
//
//
//    private static Tuple2<TestDataBreakEndVariants, TestDataBreakEndVariants> forInterChromosomeSimple_plus() {
//
//    }
//    private static Tuple2<TestDataBreakEndVariants, TestDataBreakEndVariants> forInterChromosomeSimple_minus() {
//
//    }
//
//
//    private static Tuple2<TestDataBreakEndVariants, TestDataBreakEndVariants> forInterChromosomeStrandSwitch_plus() {
//
//    }
//    private static Tuple2<TestDataBreakEndVariants, TestDataBreakEndVariants> forInterChromosomeStrandSwitch_minus() {
//
//    }
}
