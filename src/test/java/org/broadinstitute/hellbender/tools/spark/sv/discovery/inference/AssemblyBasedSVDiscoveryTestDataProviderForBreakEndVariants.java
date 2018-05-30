package org.broadinstitute.hellbender.tools.spark.sv.discovery.inference;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.TextCigarCodec;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.samtools.util.StringUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFConstants;
import org.broadinstitute.hellbender.engine.datasources.ReferenceMultiSource;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.BreakEndVariantType;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.SvType;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.TestUtilsForAssemblyBasedSVDiscovery;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AlignmentInterval;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.ContigAlignmentsModifier;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.StrandSwitch;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.testng.Assert;
import org.testng.annotations.Test;
import scala.Tuple2;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import static org.broadinstitute.hellbender.tools.spark.sv.discovery.TestUtilsForAssemblyBasedSVDiscovery.getReverseComplimentCopy;
import static org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AssemblyContigWithFineTunedAlignments.NO_GOOD_MAPPING_TO_NON_CANONICAL_CHROMOSOME;
import static org.broadinstitute.hellbender.tools.spark.sv.discovery.inference.AssemblyBasedSVDiscoveryTestDataProvider.EMPTY_BYTE_ARRAY;
import static org.broadinstitute.hellbender.tools.spark.sv.discovery.inference.AssemblyBasedSVDiscoveryTestDataProvider.makeBND;
import static org.broadinstitute.hellbender.tools.spark.sv.discovery.inference.AssemblyBasedSVDiscoveryTestDataProvider.makeBNDType;
import static org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants.*;

public class AssemblyBasedSVDiscoveryTestDataProviderForBreakEndVariants {

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
            return TestUtilsForAssemblyBasedSVDiscovery.b38_seqDict_chr20_chr21;
        }

        @Override
        public ReferenceMultiSource getAppropriateRef() {
            return TestUtilsForAssemblyBasedSVDiscovery.b38_reference_chr20_chr21;
        }

        @Override
        public Class<? extends BreakpointsInference> getAppropriateBreakpointInferencer() {
            return inferencer;
        }
    }

    public static final boolean testDataInitialized;

    public static final TestDataBreakEndVariants forIntraChromosomeRefOrderSwapWithHomology_plus;
    public static final TestDataBreakEndVariants forIntraChromosomeRefOrderSwapWithHomology_minus;
    public static final TestDataBreakEndVariants forIntraChromosomeRefOrderSwapWithInsertion_plus;
    public static final TestDataBreakEndVariants forIntraChromosomeRefOrderSwapWithInsertion_minus;
    public static final TestDataBreakEndVariants forInterChromosomeNoStrandSwitchAndUpstreamLocIsFirstInPartner_plus;
    public static final TestDataBreakEndVariants forInterChromosomeNoStrandSwitchAndUpstreamLocIsFirstInPartner_minus;
    public static final TestDataBreakEndVariants forInterChromosomeNoStrandSwitchAndUpstreamLocIsSecondInPartner_plus;
    public static final TestDataBreakEndVariants forInterChromosomeNoStrandSwitchAndUpstreamLocIsSecondInPartner_minus;
    public static final TestDataBreakEndVariants forInterChromosomeStrandSwitch55_plus;
    public static final TestDataBreakEndVariants forInterChromosomeStrandSwitch55_minus;
    public static final TestDataBreakEndVariants forInterChromosomeStrandSwitch33_plus;
    public static final TestDataBreakEndVariants forInterChromosomeStrandSwitch33_minus;

    static {
        Tuple2<TestDataBreakEndVariants, TestDataBreakEndVariants> pair = forIntraChromosomeRefOrderSwapWithHomology();
        forIntraChromosomeRefOrderSwapWithHomology_plus = pair._1;
        forIntraChromosomeRefOrderSwapWithHomology_minus = pair._2;

        pair = forIntraChromosomeRefOrderSwapWithInsertion();
        forIntraChromosomeRefOrderSwapWithInsertion_plus = pair._1;
        forIntraChromosomeRefOrderSwapWithInsertion_minus = pair._2;

        pair = forInterChromosomeNoStrandSwitchAndUpstreamLocIsFirstInPartner();
        forInterChromosomeNoStrandSwitchAndUpstreamLocIsFirstInPartner_plus = pair._1;
        forInterChromosomeNoStrandSwitchAndUpstreamLocIsFirstInPartner_minus = pair._2;

        pair = forInterChromosomeNoStrandSwitchAndUpstreamLocIsSecondInPartner();
        forInterChromosomeNoStrandSwitchAndUpstreamLocIsSecondInPartner_plus = pair._1;
        forInterChromosomeNoStrandSwitchAndUpstreamLocIsSecondInPartner_minus = pair._2;

        pair = forInterChromosomeStrandSwitch55();
        forInterChromosomeStrandSwitch55_plus = pair._1;
        forInterChromosomeStrandSwitch55_minus = pair._2;

        pair = forInterChromosomeStrandSwitch33();
        forInterChromosomeStrandSwitch33_plus = pair._1;
        forInterChromosomeStrandSwitch33_minus = pair._2;

        testDataInitialized = true;
    }

    public static List<TestDataBreakEndVariants> getAllTestData() {
        final List<TestDataBreakEndVariants> result = new ArrayList<>(40);
        for (final Tuple2<TestDataBreakEndVariants, TestDataBreakEndVariants> pair : getAllTestDataPaired()) {
            result.add(pair._1);result.add(pair._2);
        }
        return result;
    }

    public static List<Tuple2<TestDataBreakEndVariants, TestDataBreakEndVariants>> getAllTestDataPaired() {
        final List<Tuple2<TestDataBreakEndVariants, TestDataBreakEndVariants>> result = new ArrayList<>();
        result.add(forIntraChromosomeRefOrderSwapWithHomology());
        result.add(forIntraChromosomeRefOrderSwapWithInsertion());
        result.add(forInterChromosomeNoStrandSwitchAndUpstreamLocIsFirstInPartner());
        result.add(forInterChromosomeNoStrandSwitchAndUpstreamLocIsSecondInPartner());
        result.add(forInterChromosomeStrandSwitch55());
        result.add(forInterChromosomeStrandSwitch33());
        return result;
    }

    private static VariantContextBuilder addStandardAttributes(final VariantContextBuilder inputBuilder, final String contigName,
                                                               final int mapQ, final int alignLength,
                                                               final String homSeq, final String insSeq, final String mateID) {
        VariantContextBuilder builder = inputBuilder
                .attribute(TOTAL_MAPPINGS, 1)
                .attribute(CONTIG_NAMES, contigName)
                .attribute(MAPPING_QUALITIES, mapQ).attribute(HQ_MAPPINGS, mapQ >= 60 ? 1 : 0)
                .attribute(ALIGN_LENGTHS, alignLength).attribute(MAX_ALIGN_LENGTH, alignLength)
                .attribute(BND_MATEID_STR, mateID);
        if (!homSeq.isEmpty()) {
            builder.attribute(HOMOLOGY, homSeq).attribute(HOMOLOGY_LENGTH, homSeq.length());
        }
        if (!insSeq.isEmpty()) {
            builder.attribute(INSERTED_SEQUENCE, insSeq).attribute(INSERTED_SEQUENCE_LENGTH, insSeq.length());
        }
        return builder;
    }

    private static Tuple2<TestDataBreakEndVariants, TestDataBreakEndVariants> forIntraChromosomeRefOrderSwapWithHomology() {

        String contigName = "forIntraChromosomeRefOrderSwapWithHomology";
        byte[] contigSequence = "TCCACAATAGGCCTCAAATCACTCTAAATATCCACTTGCAGATTCTACAAAAAGAGTGTTTCAAAACTGCTCAATCCAAAGAAAGGTTGTACCCTGTGAGATGAATGCACGCATCACAAAGTAGTTTCTCAGAATGCTTCTGCATAGTTTTTATGTGAAGATATTTCCTTCTACACTGTAGGCCTGAAAAGGCTCCAAATATCCATTTACAGATTCTAAAAAAGAGTGTTTCAAAACTGCTATATCAATAGAAACATCCAACTCTGTGAGATGAATGCACAGATCACAAAGAAGTTTTTCAGAATGCTTCTGTGTAGTTTTTATGTGAAGATATTTGATTTTCCACAGTAGGCCCCAACGAGCTCCAAATATCCACTTGCAGATTCCACAAAAAGAGTGTTTCAAAACTGCTCAATCAACAGAGACATTCAACTCTGTGAGATGAATGCACACATCACAAAGAAGTTTCTCAAAATGCTTCTGTGTAGTTTTTGTGTGAAGATATTTCATTTTCCACAGTACGCCTCAAAGCGCTCCAAATATCCACTCTCAGATTCTGTAAAAAGAGAGATTCAAAACTGCTGAATCAAAAGATAGGTTCAACACTGTCACTTCAGTGCACAACTCACAAAGATGTTTCTCAGAATGCTTCTGTGTAGTTTTTGTGTGAAGATAATTCATTTTCCACAGTACGCCTCAAAGCGCTCCAAATATCCACTCGCAGATTCTGTAAAAAGAAAGATTCAAAACTGCTGAATCAAAAGATACGTTCAACAATGTGAGCTGAATGCACACATCACAAATAAGGTTGACAGAATGCTTCTGGGTAGTTTTTATT".getBytes();
        String homology = StringUtil.bytesToString(Arrays.copyOfRange(contigSequence, 119, 201));
        AlignmentInterval firstAlignment = new AlignmentInterval(new SimpleInterval("chr20:29189283-29189483"), 1, 201, TextCigarCodec.decode("201M637S"), true, 60, 5, 180, ContigAlignmentsModifier.AlnModType.NONE);
        AlignmentInterval secondAlignment = new AlignmentInterval(new SimpleInterval("chr20:28838813-28839531"), 120, 838, TextCigarCodec.decode("119H719M"), true, 60, 25, 594, ContigAlignmentsModifier.AlnModType.NONE);
        SimpleChimera simpleChimera = new SimpleChimera(contigName, firstAlignment, secondAlignment, StrandSwitch.NO_SWITCH, true, Collections.emptyList(), NO_GOOD_MAPPING_TO_NON_CANONICAL_CHROMOSOME);
        SimpleInterval leftBreakpoint = new SimpleInterval("chr20:28838813-28838813");
        SimpleInterval rightBreakpoint = new SimpleInterval("chr20:29189401-29189401");
        final BreakpointComplications complications = new BreakpointComplications.IntraChrRefOrderSwapBreakpointComplications(homology, "");
        NovelAdjacencyAndAltHaplotype manuallyCuratedBiPathBubble = new NovelAdjacencyAndAltHaplotype(leftBreakpoint, rightBreakpoint, StrandSwitch.NO_SWITCH, complications, TypeInferredFromSimpleChimera.INTRA_CHR_REF_ORDER_SWAP, EMPTY_BYTE_ARRAY);
        final List<SvType> types = Arrays.asList(
                makeBNDType("BND_chr20_28838813_29189401_1", Allele.create("]chr20:29189401]A"), true, BreakEndVariantType.SupportedType.INTRA_CHR_REF_ORDER_SWAP),
                makeBNDType("BND_chr20_28838813_29189401_2", Allele.create("G[chr20:28838813["), false, BreakEndVariantType.SupportedType.INTRA_CHR_REF_ORDER_SWAP)
        );

        final List<VariantContext> variants = Arrays.asList(
                addStandardAttributes(makeBND(leftBreakpoint, rightBreakpoint, Allele.create("A", true), "", "", true, false, true), contigName, 60, 119, homology, "", "BND_chr20_28838813_29189401_2").make(),
                addStandardAttributes(makeBND(leftBreakpoint, rightBreakpoint, Allele.create("G", true), "", "", false, true, false), contigName, 60, 119, homology, "", "BND_chr20_28838813_29189401_1").make()
        );

        final TestDataBreakEndVariants forIntraChromosomeRefOrderSwapWithHomology_plus =
                new TestDataBreakEndVariants(firstAlignment, secondAlignment, contigName, contigSequence, true, simpleChimera, manuallyCuratedBiPathBubble, types, variants, BreakpointsInference.IntraChrRefOrderSwapBreakpointsInference.class);


        firstAlignment = new AlignmentInterval(new SimpleInterval("chr20:28838813-28839531"), 1, 719, TextCigarCodec.decode("719M119S"), false, 60, 25, 594, ContigAlignmentsModifier.AlnModType.NONE);
        secondAlignment = new AlignmentInterval(new SimpleInterval("chr20:29189283-29189483"), 638, 838, TextCigarCodec.decode("637H201M"), false, 60, 5, 180, ContigAlignmentsModifier.AlnModType.NONE);
        simpleChimera = new SimpleChimera(contigName, firstAlignment, secondAlignment, StrandSwitch.NO_SWITCH, false, Collections.emptyList(), NO_GOOD_MAPPING_TO_NON_CANONICAL_CHROMOSOME);

        final TestDataBreakEndVariants forIntraChromosomeRefOrderSwapWithHomology_minus =
                new TestDataBreakEndVariants(firstAlignment, secondAlignment, contigName, getReverseComplimentCopy(contigSequence), false, simpleChimera, manuallyCuratedBiPathBubble, types, variants, BreakpointsInference.IntraChrRefOrderSwapBreakpointsInference.class);

        return new Tuple2<>(forIntraChromosomeRefOrderSwapWithHomology_plus, forIntraChromosomeRefOrderSwapWithHomology_minus);
    }

    private static Tuple2<TestDataBreakEndVariants, TestDataBreakEndVariants> forIntraChromosomeRefOrderSwapWithInsertion() {

        String contigName = "forIntraChromosomeRefOrderSwapWithInsertion";
        byte[] contigSequence = "TGTGGTGTGGGTGTGTGCGTGTGTGTGGACTGTGTGGTGTGGGTGTGTGCGTGTGTGTGGACTGTGTGGTGTGGGTGTGTGCGTGTGTGTGTGACCGTGTGGAGTGTGTCTGTGTGCATGTGTGGGCTGTGTGGTGCGTGTGTGCTTATGTTTGGCGTGTGTGTGTGTGTGTGTGTGGACTGTGTGGTCTGTGTGTGTGTGCGTGTGTGTGGACTCTGTGGTGTGTGCATGTGTGCATGTATGTGCGTGTGTGTGACTGTGTGGTGTGTGTCTGTGTGCACGTGTGGGCTTATATTTGGTGTGTGTGCATGTGTGGACCGTGTGGTGTGTGTGTGTGCGTGTGTGGACTGTAGTGTGTGTGCACGTGTGTGTGTGCTTGTGTGTGGAGTGTGTGTGTGTGGACTGTAGTGTGTGTGCGTGACTGTGGTGTGTGTGCATGACTGTGTGGTGTGTGTGTGCATGTGTGTGGATTGTGTGGTGTGTGTGGACTGTGGGTGTGTGGTGCGTGTGTGTGCTTGTGTGTGGTGTGTGTGCGTGTGTGGGGACTGTGTGGTGCGTGTGTGTGCTTGTGTGTGGACTGTGGATTGTGTGGTGTGTGTGTGCGCACGTGTGTGTGCGTGTCTGTGTGGTGTGTGGACTGTGTGGTGTGTGTGGACTGTGGTGTGTGTGTGCGTGACTGTGTGGTGTGTGTGTGCGCGTGTGTGTGACTGCTTGGTGTGTGTGGACTGTGGTGTGTGTGGTGTGTGTGTGCTTGTGTGTGGTGTGTG".getBytes();
        String insSeq = StringUtil.bytesToString(Arrays.copyOfRange(contigSequence, 71, 617));
        AlignmentInterval firstAlignment = new AlignmentInterval(new SimpleInterval("chr20:44405003-44405064"), 10, 71, TextCigarCodec.decode("9H62M696H"), true, 32, 3, 47, ContigAlignmentsModifier.AlnModType.NONE);
        AlignmentInterval secondAlignment = new AlignmentInterval(new SimpleInterval("chr20:44404789-44404936"), 618, 767, TextCigarCodec.decode("617S66M2I82M"), true, 60, 4, 45, ContigAlignmentsModifier.AlnModType.NONE);
        SimpleChimera simpleChimera = new SimpleChimera(contigName, firstAlignment, secondAlignment, StrandSwitch.NO_SWITCH, true, Collections.emptyList(), NO_GOOD_MAPPING_TO_NON_CANONICAL_CHROMOSOME);
        SimpleInterval leftBreakpoint = new SimpleInterval("chr20:44404789-44404789");
        SimpleInterval rightBreakpoint = new SimpleInterval("chr20:44405064-44405064");
        final BreakpointComplications complications = new BreakpointComplications.IntraChrRefOrderSwapBreakpointComplications("", insSeq);
        NovelAdjacencyAndAltHaplotype manuallyCuratedBiPathBubble = new NovelAdjacencyAndAltHaplotype(leftBreakpoint, rightBreakpoint, StrandSwitch.NO_SWITCH, complications, TypeInferredFromSimpleChimera.INTRA_CHR_REF_ORDER_SWAP, EMPTY_BYTE_ARRAY);
        final List<SvType> types = Arrays.asList(
                makeBNDType("BND_chr20_44404789_44405064_1", Allele.create("]chr20:44405064]"+insSeq+"G"), true, BreakEndVariantType.SupportedType.INTRA_CHR_REF_ORDER_SWAP),
                makeBNDType("BND_chr20_44404789_44405064_2", Allele.create("G"+insSeq+"[chr20:44404789["), false, BreakEndVariantType.SupportedType.INTRA_CHR_REF_ORDER_SWAP)
        );
        final List<VariantContext> variants = Arrays.asList(
                addStandardAttributes(makeBND(leftBreakpoint, rightBreakpoint, Allele.create("G", true), insSeq, "", true, false, true), contigName, 32, 62, "", insSeq, "BND_chr20_44404789_44405064_2").make(),
                addStandardAttributes(makeBND(leftBreakpoint, rightBreakpoint, Allele.create("G", true), insSeq, "", false, true, false), contigName, 32, 62, "", insSeq, "BND_chr20_44404789_44405064_1").make()
        );

        final TestDataBreakEndVariants forIntraChromosomeRefOrderSwapWithInsertion_plus =
                new TestDataBreakEndVariants(firstAlignment, secondAlignment, contigName, contigSequence, true, simpleChimera, manuallyCuratedBiPathBubble, types, variants, BreakpointsInference.IntraChrRefOrderSwapBreakpointsInference.class);


        firstAlignment = new AlignmentInterval(new SimpleInterval("chr20:44404789-44404936"), 1, 150, TextCigarCodec.decode("82M2I66M617S"), false, 60, 4, 45, ContigAlignmentsModifier.AlnModType.NONE);
        secondAlignment = new AlignmentInterval(new SimpleInterval("chr20:44405003-44405064"), 697, 758, TextCigarCodec.decode("696H62M9H"), false, 32, 3, 47, ContigAlignmentsModifier.AlnModType.NONE);
        simpleChimera = new SimpleChimera(contigName, firstAlignment, secondAlignment, StrandSwitch.NO_SWITCH, false, Collections.emptyList(), NO_GOOD_MAPPING_TO_NON_CANONICAL_CHROMOSOME);

        final TestDataBreakEndVariants forIntraChromosomeRefOrderSwapWithInsertion_minus =
                new TestDataBreakEndVariants(firstAlignment, secondAlignment, contigName, getReverseComplimentCopy(contigSequence), false, simpleChimera, manuallyCuratedBiPathBubble, types, variants, BreakpointsInference.IntraChrRefOrderSwapBreakpointsInference.class);

        return new Tuple2<>(forIntraChromosomeRefOrderSwapWithInsertion_plus, forIntraChromosomeRefOrderSwapWithInsertion_minus);
    }

    private static Tuple2<TestDataBreakEndVariants, TestDataBreakEndVariants> forInterChromosomeNoStrandSwitchAndUpstreamLocIsFirstInPartner() {

        String contigName = "forInterChromosomeNoStrandSwitchAndUpstreamLocIsFirstInPartner";
        byte[] contigSequence = "GATTCTACAAAAAGACTGTTTCCAAACTGCTCAATCAAAAGAAAGGTTCAAACCTGTGAGATGAAAGCACACATCACTAAGAAGTTTCTCAGAATGCTTCTGTCTAGTTTTTAAGTGAAGATATTTCCTCTTTCACCATAGGTCTCAATGGGCTCAGAAATATCCCTTTGCAGATCCTACAAAATGACTGTTTCCAAACTGCTCAATCAAAACAAAGTTTCAACTCTGTCAGATGAATGCACACATCACAAAGAAATTTCTCAGAATGCTTCTGTCTAGTTTTTATGTGAAGATATTTCCTTTTCCACCATAGGCCTCAAAGCGCTCCAAATATCCACTTGAACACTCTACGAAAGAGTTTTTCAAAACTGCTCAATCAAAAGAAAGCTTCAGCTCTGTGACATAAATGCACACAGCACAAAGAAGTTTCTCAGAATGCTTGTGTCTAGTTTTTATGTGAAGATGTTTCCTTTTCCACCAGAGGCCTCAAAGCACTTCAAATATACACTTGCAGATACTGCAAAAAGAGGGTTTCAAAACTGCTCAATCAAAAGAAAGTTGGAAGTCTGTGAGATGAATGCACACATCACAAAGAAGTTTCTAAGAATGCTTCCATCTGAATTTTATGTGAGGATATTTCCTTTTTCACCATAG".getBytes();
        String homology = StringUtil.bytesToString(Arrays.copyOfRange(contigSequence, 166, 322));
        AlignmentInterval firstAlignment = new AlignmentInterval(new SimpleInterval("chr20:28766544-28766865"), 1, 322, TextCigarCodec.decode("322M332H"), true, 58, 9, 277, ContigAlignmentsModifier.AlnModType.NONE);
        AlignmentInterval secondAlignment = new AlignmentInterval(new SimpleInterval("chr21:7969682-7970168"), 167, 654, TextCigarCodec.decode("166S211M1I276M"), true, 58, 33, 310, ContigAlignmentsModifier.AlnModType.NONE);
        SimpleChimera simpleChimera = new SimpleChimera(contigName, firstAlignment, secondAlignment, StrandSwitch.NO_SWITCH, true, Collections.emptyList(), NO_GOOD_MAPPING_TO_NON_CANONICAL_CHROMOSOME);
        SimpleInterval leftBreakpoint = new SimpleInterval("chr20:28766709-28766709");
        SimpleInterval rightBreakpoint = new SimpleInterval("chr21:7969682-7969682");
        final BreakpointComplications complications = new BreakpointComplications.InterChromosomeBreakpointComplications(homology, "");
        NovelAdjacencyAndAltHaplotype manuallyCuratedBiPathBubble = new NovelAdjacencyAndAltHaplotype(leftBreakpoint, rightBreakpoint, StrandSwitch.NO_SWITCH, complications, TypeInferredFromSimpleChimera.INTER_CHR_NO_SS_WITH_LEFT_MATE_FIRST_IN_PARTNER, EMPTY_BYTE_ARRAY);
        List<SvType> types = Arrays.asList(
                makeBNDType("BND_chr20_28766709_chr21_7969682_1", Allele.create("C[chr21:7969682["), true, BreakEndVariantType.SupportedType.INTER_CHR_NO_SS_WITH_LEFT_MATE_FIRST_IN_PARTNER),
                makeBNDType("BND_chr20_28766709_chr21_7969682_2", Allele.create("]chr20:28766709]T"), false, BreakEndVariantType.SupportedType.INTER_CHR_NO_SS_WITH_LEFT_MATE_FIRST_IN_PARTNER)
        );
        List<VariantContext> variants = Arrays.asList(
                addStandardAttributes(makeBND(leftBreakpoint, rightBreakpoint, Allele.create("C", true), "", "", true, true, false), contigName, 58, 166, homology, "", "BND_chr20_28766709_chr21_7969682_2").make(),
                addStandardAttributes(makeBND(leftBreakpoint, rightBreakpoint, Allele.create("T", true), "", "", false, false, true), contigName, 58, 166, homology, "", "BND_chr20_28766709_chr21_7969682_1").make()
        );

        final TestDataBreakEndVariants forInterChromosomeNoStrandSwitchAndUpstreamLocIsFirstInPartner_plus =
                new TestDataBreakEndVariants(firstAlignment, secondAlignment, contigName, contigSequence, false, simpleChimera, manuallyCuratedBiPathBubble, types, variants, BreakpointsInference.InterChromosomeBreakpointsInference.class);


        contigSequence = "AACTTTCTTTTGATTGAGCAGTTTTGAAACCCTCTTTTTGCAGTATCTGCAAGTGTATATTTGAAGTGCTTTGAGGCCTCTGGTGGAAAAGGAAACATCTTCACATAAAAACTAGACACAAGCATTCTGAGAAACTTCTTTGTGCTGTGTGCATTTATGTCACAGAGCTGAAGCTTTCTTTTGATTGAGCAGTTTTGAAAAACTCTTTCGTAGAGTGTTCAAGTGGATATTTGGAGCGCTTTGAGGCCTATGGTGGAAAAGGAAATATCTTCACATAAAAACTAGACAGAAGCATTCTGAGAAATTTCTTTGTGATGTGTGCATTCATCTGACAGAGTTGAAACTTTGTTTTGATTGAGCAGTTTGGAAACAGTCATTTTGTAGGATCTGCAAAGGGATATTTCTGAGCCCATTGAGACCTATGGTGAAAGAGGAAATATCTTCACTTAAAAACTAGACAGAAGCATTCTGAGAAACTTCTTAGTGATGTGTGCTTTCATCTCACAGGTTTGAACCTTTCTTTTAATTGAGCAGTTTGGAAACAGTCTTTTTGTGGAATCTGCAAAGGATAATTCAAGCACTTTGAGGCCTATGGTGAAAAAGGACATATCTTCACATGAAATCTAAACAGAAGCTTTCTGAGAAACTTCTTTTTGATGAGTGCATACATCTCACAGAGGTGAAACTTTCTTTTCATTGAACAGCTTGGAAACAGTCTTTTTGTACAATCTGCAAAGGAATATTTCTGCGAAGTTAGAGGCCTATGGTGAAAAAGAAATATCTTCAGATAAAATGTAGACAGAAGTATTCTGAGAAAATTTTTTGTGATGTATCCATTCATCTCACAGAGTTGAACTTTTCTTTTGATGGAGCAGTCTGGAAACAGTCTTTTTGTAGTATCTGAAGAGGTATATGTGAG".getBytes();
        firstAlignment = new AlignmentInterval(new SimpleInterval("chr21:7969682-7970074"), 1, 394, TextCigarCodec.decode("182M1I211M525H"), false, 58, 23, 266, ContigAlignmentsModifier.AlnModType.NONE);
        secondAlignment = new AlignmentInterval(new SimpleInterval("chr20:28766185-28766865"), 239, 919, TextCigarCodec.decode("238S681M"), false, 58, 7, 646, ContigAlignmentsModifier.AlnModType.NONE);
        simpleChimera = new SimpleChimera(contigName, firstAlignment, secondAlignment, StrandSwitch.NO_SWITCH, false, Collections.emptyList(), NO_GOOD_MAPPING_TO_NON_CANONICAL_CHROMOSOME);
        variants = Arrays.asList(
                addStandardAttributes(makeBND(leftBreakpoint, rightBreakpoint, Allele.create("C", true), "", "", true, true, false), contigName, 58, 237, homology, "", "BND_chr20_28766709_chr21_7969682_2").make(),
                addStandardAttributes(makeBND(leftBreakpoint, rightBreakpoint, Allele.create("T", true), "", "", false, false, true), contigName, 58, 237, homology, "", "BND_chr20_28766709_chr21_7969682_1").make()
        );
        final TestDataBreakEndVariants forInterChromosomeNoStrandSwitchAndUpstreamLocIsFirstInPartner_minus =
                new TestDataBreakEndVariants(firstAlignment, secondAlignment, contigName, contigSequence, true, simpleChimera, manuallyCuratedBiPathBubble, types, variants, BreakpointsInference.InterChromosomeBreakpointsInference.class);

        return new Tuple2<>(forInterChromosomeNoStrandSwitchAndUpstreamLocIsFirstInPartner_plus, forInterChromosomeNoStrandSwitchAndUpstreamLocIsFirstInPartner_minus);
    }

    private static Tuple2<TestDataBreakEndVariants, TestDataBreakEndVariants> forInterChromosomeNoStrandSwitchAndUpstreamLocIsSecondInPartner() {

        String contigName = "forInterChromosomeNoStrandSwitchAndUpstreamLocIsSecondInPartner";
        byte[] contigSequence = "CAAAGAAGTTTCTCATAATGCTTCTGTCTACTTGTTGTGTGAAGATACTTCCTATTTCACCATAGGCCTCAAAGGCCTCACAAATATCCCTTTGCAGACTCTACAAAAAGACAGTTTTTGAACTACTCCATGAAAAGAAAGTTTCAACTCTGAGTCGAATTTCACACAGAATAATTTAGTTTCTCATAATGTTTCTGTCTAGTTTTTATGTGACGGTATTTCCTTTTTCACCATAGGCCTCAAACCGCTCACAAATATCCCTCTGCAGATACAACAAAAGGACAGTTTCCAAACTGCTAAATCAAAAGATATGCTCAACAACGTGAAATGAATGCACATGTCACAAAGAAGTTTCTCAGAATGCTTCTGTCTAGTTTTAATGTGAAGATATTTACTTTTTCTCCATAGGCCCCAAAGCACTCCAAACGTCCATTTGCAGATTCCACAAAAAGACTGTTTCCAAACTTCTCAATCAAAAGAAAGGTTCAATTCTGTGAGATGAAAGCACACATCACAAAGAAGTTTCTCAGAAAGCTTCTGTCTAGTTTTTATGTGAGGATATTTTATATTTCACCACAGGCCTCAATGGACTCAAAAATATCCCTATCAGATTCTACAAAAAGACTGTGTCCGAACTTCTCAATGAAGAGAAACTTTCAACTCTGTGAGGTGAATGCACACATAAAAAAGAAGTTTCTCAGAATGCTTCTGTCTAGTTTTTATGTGAAGATACTACTTTTTCACCATAGACCTCAAACCGCTCAGAAATATCCCTTTGCAGATTGTACAAAAAGACTGTTTCCAAACTGCTCAATAAAAAGAAAGATTCATCTCTGTGAGATAAATGCAAAAATAATAAAGAAGTTTCTCAGAATGCTTCTGTCTAGTTTTTATGTGAAGATATTTCCTTTTTCACCATAGTCCTTAAACTGCTCACAAATATCTCTCTGCAGATACTACAAAAAGACTGTTTCCAAACTGCTCCATGAAAAGAAAGGTTCAACTCTGTGAGACGAATGCACACATCACAAAGAATATTCTCAGAATGATTCCATCTAAT".getBytes();
        String homology = "";
        String insSeq = StringUtil.bytesToString(Arrays.copyOfRange(contigSequence, 152, 192));
        AlignmentInterval firstAlignment = new AlignmentInterval(new SimpleInterval("chr21:10806570-10806721"), 1, 152, TextCigarCodec.decode("152M908H"), true, 22, 14, 82, ContigAlignmentsModifier.AlnModType.NONE);
        AlignmentInterval secondAlignment = new AlignmentInterval(new SimpleInterval("chr20:29286550-29287416"), 193, 1060, TextCigarCodec.decode("192S604M1I263M"), true, 60, 33, 690, ContigAlignmentsModifier.AlnModType.NONE);
        SimpleChimera simpleChimera = new SimpleChimera(contigName, firstAlignment, secondAlignment, StrandSwitch.NO_SWITCH, true, Collections.emptyList(), NO_GOOD_MAPPING_TO_NON_CANONICAL_CHROMOSOME);
        SimpleInterval leftBreakpoint = new SimpleInterval("chr20:29286550-29286550");
        SimpleInterval rightBreakpoint = new SimpleInterval("chr21:10806721-10806721");
        final BreakpointComplications complications = new BreakpointComplications.InterChromosomeBreakpointComplications(homology, insSeq);
        NovelAdjacencyAndAltHaplotype manuallyCuratedBiPathBubble = new NovelAdjacencyAndAltHaplotype(leftBreakpoint, rightBreakpoint, StrandSwitch.NO_SWITCH, complications, TypeInferredFromSimpleChimera.INTER_CHR_NO_SS_WITH_LEFT_MATE_SECOND_IN_PARTNER, EMPTY_BYTE_ARRAY);
        List<SvType> types = Arrays.asList(
                makeBNDType("BND_chr20_29286550_chr21_10806721_1", Allele.create("]chr21:10806721]"+insSeq+"T"), true, BreakEndVariantType.SupportedType.INTER_CHR_NO_SS_WITH_LEFT_MATE_SECOND_IN_PARTNER),
                makeBNDType("BND_chr20_29286550_chr21_10806721_2", Allele.create("G"+insSeq+"[chr20:29286550["), false, BreakEndVariantType.SupportedType.INTER_CHR_NO_SS_WITH_LEFT_MATE_SECOND_IN_PARTNER)
        );
        List<VariantContext> variants = Arrays.asList(
                addStandardAttributes(makeBND(leftBreakpoint, rightBreakpoint, Allele.create("T", true), insSeq, "", true, false, true), contigName, 22, 152, homology, insSeq, "BND_chr20_29286550_chr21_10806721_2").make(),
                addStandardAttributes(makeBND(leftBreakpoint, rightBreakpoint, Allele.create("G", true), insSeq, "", false, true, false), contigName, 22, 152, homology, insSeq, "BND_chr20_29286550_chr21_10806721_1").make()
        );

        final TestDataBreakEndVariants forInterChromosomeNoStrandSwitchAndUpstreamLocIsSecondInPartner_plus =
                new TestDataBreakEndVariants(firstAlignment, secondAlignment, contigName, contigSequence, true, simpleChimera, manuallyCuratedBiPathBubble, types, variants, BreakpointsInference.InterChromosomeBreakpointsInference.class);


        firstAlignment = new AlignmentInterval(new SimpleInterval("chr20:29286550-29287416"), 1, 868, TextCigarCodec.decode("192S604M1I263M"), false, 60, 33, 690, ContigAlignmentsModifier.AlnModType.NONE);
        secondAlignment = new AlignmentInterval(new SimpleInterval("chr21:10806570-10806721"), 909, 1060, TextCigarCodec.decode("152M908H"), false, 22, 14, 82, ContigAlignmentsModifier.AlnModType.NONE);
        simpleChimera = new SimpleChimera(contigName, firstAlignment, secondAlignment, StrandSwitch.NO_SWITCH, false, Collections.emptyList(), NO_GOOD_MAPPING_TO_NON_CANONICAL_CHROMOSOME);
        final TestDataBreakEndVariants forInterChromosomeNoStrandSwitchAndUpstreamLocIsSecondInPartner_minus =
                new TestDataBreakEndVariants(firstAlignment, secondAlignment, contigName, getReverseComplimentCopy(contigSequence), false, simpleChimera, manuallyCuratedBiPathBubble, types, variants, BreakpointsInference.InterChromosomeBreakpointsInference.class);

        return new Tuple2<>(forInterChromosomeNoStrandSwitchAndUpstreamLocIsSecondInPartner_plus, forInterChromosomeNoStrandSwitchAndUpstreamLocIsSecondInPartner_minus);
    }

    private static Tuple2<TestDataBreakEndVariants, TestDataBreakEndVariants> forInterChromosomeStrandSwitch55() {

        String contigName = "forInterChromosomeStrandSwitch55";
        byte[] contigSequence = "ATTCTGAGAAACTTCATTTTGATGTGTGCATTCATCTTCCAGAGTTGAAACTTTCTTTTGATTGTGTAGTTTTGAAACACTCTTTTTGTAGAATCTGCAAGGGGGTATTTGTAGGGATTTGAAGCCTATTGTTGAAAAGGTAATATCTTCACATAAAAACTACATAGGATCATTCGGAGAAACTTCTTTGTGATGTGTGCATTCAACTCACAGAGTTGAACCTATCTTTTTTTTTTTAATGTCTGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTGGAGCAATTTGAGGCCTAAGGTGGAAAAGGAAATATTTTCACATAAAAACTAGACAGAAGAATTCTGTGAAACTTGTTCAGGACCTGTGCATTCATCTTACAGATTTGAATCTTTCTTTTGATTGAGCAGTTTGGAAACACTGTTTTTGTAGAATCTTCAGGTGGACATTCAGAGCACTTTGTGTCCTATGGTAGAAAAGGAAATATCTTCATA".getBytes();
        String homology = "";
        String insSeq = StringUtil.bytesToString(getReverseComplimentCopy(Arrays.copyOfRange(contigSequence, 230, 292)));
        AlignmentInterval firstAlignment = new AlignmentInterval(new SimpleInterval("chr21:10784600-10784829"), 1, 230, TextCigarCodec.decode("230M279S"), true, 60, 12, 170, ContigAlignmentsModifier.AlnModType.NONE);
        AlignmentInterval secondAlignment = new AlignmentInterval(new SimpleInterval("chr20:28817762-28817977"), 293, 509, TextCigarCodec.decode("292H93M1I123M"), false, 60, 11, 149, ContigAlignmentsModifier.AlnModType.NONE);
        SimpleChimera simpleChimera = new SimpleChimera(contigName, firstAlignment, secondAlignment, StrandSwitch.FORWARD_TO_REVERSE, false, Collections.emptyList(), NO_GOOD_MAPPING_TO_NON_CANONICAL_CHROMOSOME);
        SimpleInterval leftBreakpoint = new SimpleInterval("chr20:28817977-28817977");
        SimpleInterval rightBreakpoint = new SimpleInterval("chr21:10784829-10784829");
        final BreakpointComplications complications = new BreakpointComplications.InterChromosomeBreakpointComplications(homology, insSeq);
        NovelAdjacencyAndAltHaplotype manuallyCuratedBiPathBubble = new NovelAdjacencyAndAltHaplotype(leftBreakpoint, rightBreakpoint, StrandSwitch.FORWARD_TO_REVERSE, complications, TypeInferredFromSimpleChimera.INTER_CHR_STRAND_SWITCH_55, EMPTY_BYTE_ARRAY);
        final List<SvType> types = Arrays.asList(
                makeBNDType("BND_chr20_28817977_chr21_10784829_1", Allele.create("A"+insSeq+"]chr21:10784829]"), true, BreakEndVariantType.SupportedType.INTER_CHR_STRAND_SWITCH_55),
                makeBNDType("BND_chr20_28817977_chr21_10784829_2", Allele.create("T"+ SequenceUtil.reverseComplement(insSeq) +"]chr20:28817977]"), false, BreakEndVariantType.SupportedType.INTER_CHR_STRAND_SWITCH_55)
        );
        final List<VariantContext> variants = Arrays.asList(
                addStandardAttributes(makeBND(leftBreakpoint, rightBreakpoint, Allele.create("A", true), insSeq, "", true, true, true), contigName, 60, 216, homology, insSeq, "BND_chr20_28817977_chr21_10784829_2").make(),
                addStandardAttributes(makeBND(leftBreakpoint, rightBreakpoint, Allele.create("T", true), SequenceUtil.reverseComplement(insSeq), "", false, true, true), contigName, 60, 216, homology, insSeq, "BND_chr20_28817977_chr21_10784829_1").make()
        );

        final TestDataBreakEndVariants forInterChromosomeStrandSwitch55_plus =
                new TestDataBreakEndVariants(firstAlignment, secondAlignment, contigName, contigSequence, true, simpleChimera, manuallyCuratedBiPathBubble, types, variants, BreakpointsInference.InterChromosomeBreakpointsInference.class);


        firstAlignment = new AlignmentInterval(new SimpleInterval("chr20:28817762-28817977"), 1, 217, TextCigarCodec.decode("123M1I93M292H"), true, 60, 11, 149, ContigAlignmentsModifier.AlnModType.NONE);
        secondAlignment = new AlignmentInterval(new SimpleInterval("chr21:10784600-10784829"), 280, 509, TextCigarCodec.decode("279S230M"), false, 60, 12, 170, ContigAlignmentsModifier.AlnModType.NONE);
        simpleChimera = new SimpleChimera(contigName, firstAlignment, secondAlignment, StrandSwitch.FORWARD_TO_REVERSE, true, Collections.emptyList(), NO_GOOD_MAPPING_TO_NON_CANONICAL_CHROMOSOME);

        final TestDataBreakEndVariants forInterChromosomeStrandSwitch55_minus =
                new TestDataBreakEndVariants(firstAlignment, secondAlignment, contigName, getReverseComplimentCopy(contigSequence), false, simpleChimera, manuallyCuratedBiPathBubble, types, variants, BreakpointsInference.InterChromosomeBreakpointsInference.class);

        return new Tuple2<>(forInterChromosomeStrandSwitch55_plus, forInterChromosomeStrandSwitch55_minus);
    }

    private static Tuple2<TestDataBreakEndVariants, TestDataBreakEndVariants> forInterChromosomeStrandSwitch33() {

        String contigName = "forInterChromosomeStrandSwitch33";
        byte[] contigSequence = "AACTGCACAATGAAAAGTAAGTTTCAACTCTGTGAGATGAAAGCACACATCATGAAGAAGTTTGTCAGAATGCTTCTGTCTACTTTGTATTTGAAGATATTTTTTCTTTTCCACTTTAGACCTAAAAGCGCTCCAAATGTCCACTTGCAGATTCTACAAAAAGAGAGTTTCAAAGCTGCTCAATGAAAATAAAGTTTCAACTCTGTGAGATGAATGCACACATCACAAAAAAGTTTGTCAGAATGCATCTGTCTAGTTTTTATGTGAAGATACTTCCTTTTCCACCACAGCCCCCAAAGCACTCCAAATATCCACTTGCAGGTTCTACAAAAAGAGTGTTTGCAAACTGCTCAATGAAAAGTAAGGTTCAAATTTGTGAGATGAATGCACACATCACAAAGAAGTTTGTCAGAATGCTTTGGTCTAGTTTTTATGTGAAGATATTTCCTTTTCCACCATAGGCCTCAAAGTGCGAAAAATGTCCATTGCAGATTCTACAGAAAGAGTGTTTCCAATCTGTTCAATGAAATGGAGGGTTCAACTCTGTGAGATGAATGCACACATCAAAAAGAAGTTTGTCAGAATGCTTCTGTCCCATTATTACATGAAAATATTTCGTTGTCCACCATAGGCCTCAAAGCTCTCCAAATGTCAACTTGCAGATTCTACAAAAAGCGTGTTTCAAAGCTGCTCAATCAAAAGGAAGTTTTAACTCTGTGAGATGATTGCACATGTCACAAAGAAGTTTGTCAGAACGCTTCTGTTTAGTTTTTATGTGAAGATATTTCCT".getBytes();
        String homology = SequenceUtil.reverseComplement(StringUtil.bytesToString(Arrays.copyOfRange(contigSequence, 375, 530)));
        String insSeq = "";
        AlignmentInterval firstAlignment = new AlignmentInterval(new SimpleInterval("chr20:29204769-29205298"), 1, 530, TextCigarCodec.decode("530M260S"), false, 58, 2, 520, ContigAlignmentsModifier.AlnModType.NONE);
        AlignmentInterval secondAlignment = new AlignmentInterval(new SimpleInterval("chr21:7229293-7229707"), 376, 790, TextCigarCodec.decode("375H110M1D44M1I260M"), true, 25, 57, 105, ContigAlignmentsModifier.AlnModType.NONE);
        SimpleChimera simpleChimera = new SimpleChimera(contigName, firstAlignment, secondAlignment, StrandSwitch.REVERSE_TO_FORWARD, true, Collections.emptyList(), NO_GOOD_MAPPING_TO_NON_CANONICAL_CHROMOSOME);
        SimpleInterval leftBreakpoint = new SimpleInterval("chr20:29204769-29204769");
        SimpleInterval rightBreakpoint = new SimpleInterval("chr21:7229448-7229448");
        final BreakpointComplications complications = new BreakpointComplications.InterChromosomeBreakpointComplications(homology, insSeq);
        NovelAdjacencyAndAltHaplotype manuallyCuratedBiPathBubble = new NovelAdjacencyAndAltHaplotype(leftBreakpoint, rightBreakpoint, StrandSwitch.REVERSE_TO_FORWARD, complications, TypeInferredFromSimpleChimera.INTER_CHR_STRAND_SWITCH_33, EMPTY_BYTE_ARRAY);
        final List<SvType> types = Arrays.asList(
                makeBNDType("BND_chr20_29204769_chr21_7229448_1", Allele.create("[chr21:7229448[A"), true, BreakEndVariantType.SupportedType.INTER_CHR_STRAND_SWITCH_33),
                makeBNDType("BND_chr20_29204769_chr21_7229448_2", Allele.create("[chr20:29204769[A"), false, BreakEndVariantType.SupportedType.INTER_CHR_STRAND_SWITCH_33)
        );
        final List<VariantContext> variants = Arrays.asList(
                addStandardAttributes(makeBND(leftBreakpoint, rightBreakpoint, Allele.create("A", true), insSeq, "", true, false, false), contigName, 25, 260, homology, insSeq, "BND_chr20_29204769_chr21_7229448_2").make(),
                addStandardAttributes(makeBND(leftBreakpoint, rightBreakpoint, Allele.create("A", true), insSeq, "", false, false, false), contigName, 25, 260, homology, insSeq, "BND_chr20_29204769_chr21_7229448_1").make()
        );

        final TestDataBreakEndVariants forInterChromosomeStrandSwitch33_plus =
                new TestDataBreakEndVariants(firstAlignment, secondAlignment, contigName, contigSequence, false, simpleChimera, manuallyCuratedBiPathBubble, types, variants, BreakpointsInference.InterChromosomeBreakpointsInference.class);


        firstAlignment = new AlignmentInterval(new SimpleInterval("chr21:7229293-7229707"), 1, 415, TextCigarCodec.decode("260M1I44M1D110M375H"), false, 25, 57, 105, ContigAlignmentsModifier.AlnModType.NONE);
        secondAlignment = new AlignmentInterval(new SimpleInterval("chr20:29204769-29205298"), 261, 790, TextCigarCodec.decode("260S530M"), true, 58, 2, 520, ContigAlignmentsModifier.AlnModType.NONE);
        simpleChimera = new SimpleChimera(contigName, firstAlignment, secondAlignment, StrandSwitch.REVERSE_TO_FORWARD, false, Collections.emptyList(), NO_GOOD_MAPPING_TO_NON_CANONICAL_CHROMOSOME);

        final TestDataBreakEndVariants forInterChromosomeStrandSwitch33_minus =
                new TestDataBreakEndVariants(firstAlignment, secondAlignment, contigName, getReverseComplimentCopy(contigSequence), true, simpleChimera, manuallyCuratedBiPathBubble, types, variants, BreakpointsInference.InterChromosomeBreakpointsInference.class);

        return new Tuple2<>(forInterChromosomeStrandSwitch33_plus, forInterChromosomeStrandSwitch33_minus);
    }
}
