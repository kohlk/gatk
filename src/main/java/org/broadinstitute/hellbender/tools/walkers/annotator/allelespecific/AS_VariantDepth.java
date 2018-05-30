package org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific;



        import htsjdk.variant.variantcontext.Allele;
        import htsjdk.variant.variantcontext.Genotype;
        import htsjdk.variant.variantcontext.GenotypesContext;
        import htsjdk.variant.variantcontext.VariantContext;
        import htsjdk.variant.vcf.VCFInfoHeaderLine;
        import org.broadinstitute.hellbender.tools.walkers.annotator.InfoFieldAnnotation;
        import org.broadinstitute.gatk.utils.MathUtils;
        import org.broadinstitute.hellbender.engine.AlignmentContext;
        import org.broadinstitute.hellbender.engine.ReferenceContext;
        import org.broadinstitute.hellbender.utils.genotyper.ReadLikelihoods;
        import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;
        import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
        import org.broadinstitute.hellbender.utils.variant.GATKVCFHeaderLines;

        import java.util.*;

/**
 * Created by gauthier on 8/7/17.
 */
public class AS_VariantDepth extends InfoFieldAnnotation implements AS_StandardAnnotation, ReducibleAnnotation {
    protected final String printDelim = ",";
    protected final String printFormat = "%d";
    private final String splitDelim = ",";

    @Override
    public List<String> getKeyNames() { return Arrays.asList(GATKVCFConstants.AS_VARIANT_DEPTH_KEY); }

    @Override
    public String getRawKeyName() { return GATKVCFConstants.AS_VARIANT_DEPTH_KEY; }

    public List<VCFInfoHeaderLine> getDescriptions() {
        return Arrays.asList(GATKVCFHeaderLines.getInfoLine(getKeyNames().get(0)));
    }

    public Map<String, Object> annotate(final RefMetaDataTracker tracker,
                                        final AnnotatorCompatible walker,
                                        final ReferenceContext ref,
                                        final Map<String, AlignmentContext> stratifiedContexts,
                                        final VariantContext vc,
                                        final Map<String, PerReadAlleleLikelihoodMap> perReadAlleleLikelihoodMap ) {
        return null;
    }

    @Override
    public Map<String, Object> annotateRawData(RefMetaDataTracker tracker,
                                               AnnotatorCompatible walker,
                                               ReferenceContext ref, Map<String, AlignmentContext> stratifiedContexts,
                                               VariantContext vc,
                                               Map<String, PerReadAlleleLikelihoodMap> stratifiedPerReadAlleleLikelihoodMap){
        //TODO: is used by AS_QualByDepth, but there's no real dependency at this level, just an error later on

        final GenotypesContext genotypes = vc.getGenotypes();
        if ( genotypes == null || genotypes.isEmpty() )
            return null;

        final Map<String, Object> annotations = new HashMap<>();
        final AlleleSpecificAnnotationData<Number> myRawData = new AlleleSpecificAnnotationData<>(vc.getAlleles(), null);
        calculateRawData(vc, stratifiedPerReadAlleleLikelihoodMap, myRawData);
        String annotationString = makeRawAnnotationString(vc.getAlleles(), myRawData.getAttributeMap());
        annotations.put(getRawKeyName(), annotationString);
        return annotations;
    }

    @Override
    public Map<String, Object> combineRawData(List<Allele> allelesList, List<? extends ReducibleAnnotationData> listOfRawData) {
        //VC already contains merged alleles from ReferenceConfidenceVariantContextMerger
        ReducibleAnnotationData combinedData = new AlleleSpecificAnnotationData(allelesList, null);

        for (final ReducibleAnnotationData currentValue : listOfRawData) {
            parseRawDataString(currentValue);
            combineAttributeMap(currentValue, combinedData);

        }
        final Map<String, Object> annotations = new HashMap<>();
        String annotationString = makeRawAnnotationString(allelesList, combinedData.getAttributeMap());
        annotations.put(getRawKeyName(), annotationString);
        return annotations;
    }

    protected void parseRawDataString(final ReducibleAnnotationData<Number> myData) {
        final String rawDataString = myData.getRawData();
        //get per-allele data by splitting on allele delimiter
        final String[] rawDataPerAllele = rawDataString.split(splitDelim);
        for (int i=0; i<rawDataPerAllele.length; i++) {
            final String alleleData = rawDataPerAllele[i];
            myData.putAttribute(myData.getAlleles().get(i), Integer.parseInt(alleleData));
        }
    }

    public void combineAttributeMap(final ReducibleAnnotationData<Number> toAdd, final ReducibleAnnotationData<Number> combined) {
        //check that alleles match
        for (final Allele currentAllele : combined.getAlleles()){
            //combined is initialized with all alleles, but toAdd might have only a subset
            if(toAdd.getAttribute(currentAllele) == null)
                continue;
            if (toAdd.getAttribute(currentAllele) != null && combined.getAttribute(currentAllele) != null) {
                combined.putAttribute(currentAllele, (int) combined.getAttribute(currentAllele) + (int) toAdd.getAttribute(currentAllele));
            }
            else
                combined.putAttribute(currentAllele, toAdd.getAttribute(currentAllele));
        }
    }

    @Override
    public Map<String, Object> finalizeRawData(VariantContext vc, VariantContext originalVC) {
        //TODO: for now, do all the calcs in finalize
        final GenotypesContext genotypes = vc.getGenotypes();
        if ( genotypes == null || genotypes.isEmpty() )
            return null;

        final List<Integer> standardDepth = getAlleleDepths(genotypes);
        if (standardDepth == null) //all no-calls and homRefs
            return null;

        final Map<String, Object> map = new HashMap<>();
        map.put(getKeyNames().get(0), AnnotationUtils.encodeIntValueList(standardDepth));
        return map;

        /*
        //we need to use the AS_QUAL value that was added to the VC by the GenotypingEngine
        if ( !vc.hasAttribute(getRawKeyName()) )
            return null;

        //TODO: we can take this out
        //Parse the VC's allele-specific qual values
        List<Object> alleleVarDPObjList = vc.getAttributeAsList(getRawKeyName());
        if (alleleVarDPObjList.size() != vc.getNAlleles() -1)
            throw new IllegalStateException("Number of " + getRawKeyName() +" values doesn't match the number of alternate alleles.");
        List<Integer> alleleVarDPList = new ArrayList<>();
        for (final Object obj : alleleVarDPObjList) {
            alleleVarDPList.add(Integer.parseInt(obj.toString()));
        }

        final Map<String, Object> map = new HashMap<>();
        map.put(getKeyNames().get(0), AnnotationUtils.encodeIntValueList(alleleVarDPList));
        return map;*/
    }

    public void calculateRawData(final VariantContext vc,
                                 final Map<String, PerReadAlleleLikelihoodMap> stratifiedPerReadAlleleLikelihoodMap,
                                 final ReducibleAnnotationData myRawData) {
        List<Integer> alleleDepths = getAlleleDepths(vc.getGenotypes());

        for(int i = 0; i < vc.getNAlleles(); i++) {
            Allele current = vc.getAlleles().get(i);
            Integer dp = alleleDepths.get(i);
            myRawData.putAttribute(current, dp);
        }
    }

    public String makeRawAnnotationString(final List<Allele> vcAlleles, final Map<Allele, Number> perAlleleValues) {
        String annotationString = "";
        for (final Allele current : vcAlleles) {
            if(current.isReference())
                continue;
            if (!annotationString.isEmpty())
                annotationString += printDelim;
            if(perAlleleValues.get(current) != null)
                annotationString += String.format(printFormat,perAlleleValues.get(current));
            else
                annotationString += String.format(printFormat, 0);
        }
        return annotationString;
    }

    //returns a list with the same dimension as the AD vector, which should be all alleles including reference
    private List<Integer> getAlleleDepths(final GenotypesContext genotypes) {
        int numAlleles = -1;
        for (final Genotype genotype : genotypes) {
            if (genotype.hasAD()) {
                numAlleles = genotype.getAD().length;
                break;
            }
        }
        if (numAlleles == -1) //no genotypes have AD
            return null;
        Integer[] alleleDepths = new Integer[numAlleles];
        for (int i = 0; i < alleleDepths.length; i++) {
            alleleDepths[i] = 0;
        }
        for (final Genotype genotype : genotypes) {
            // we care only about genotypes with variant alleles
            if ( !genotype.isHet() && !genotype.isHomVar() )
                continue;

            // if we have the AD values for this sample, let's make sure that the variant depth is greater than 1!
            if ( genotype.hasAD() ) {
                final int[] AD = genotype.getAD();
                final int totalADdepth = (int) MathUtils.sum(AD);
                if ( totalADdepth - AD[0] > 1 ) {
                    for (int i = 0; i < AD.length; i++) {
                        alleleDepths[i] += AD[i];
                    }
                }
            }
        }
        return Arrays.asList(alleleDepths);
    }
}

