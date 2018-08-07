/*
 * CenancestorTreeLikelihood.java
 *
 * Modified from TreeLikelihood.java by DM
 * 
 * Copyright (c) 2002-2015 Alexei Drummond, Andrew Rambaut and Marc Suchard
 *
 * This file is part of BEAST.
 * See the NOTICE file distributed with this work for additional
 * information regarding copyright ownership and licensing.
 *
 * BEAST is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 *  BEAST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package PISCA;

import dr.evolution.alignment.AscertainedSitePatterns;
import dr.evolution.alignment.PatternList;
import dr.evolution.alignment.SitePatterns;
import dr.evolution.datatype.DataType;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.util.Taxon;
import dr.evolution.util.TaxonList;
import dr.evomodel.sitemodel.SiteModel;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.treelikelihood.AbstractTreeLikelihood;
import dr.evomodel.treelikelihood.TipStatesModel;
import dr.evomodel.treelikelihood.LikelihoodCore;
import dr.inference.model.Model;
import dr.inference.model.Statistic;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;

import java.util.logging.Logger;

/**
 * CenancestorTreeLikelihoodModel - implements a Likelihood Function for sequences on a tree with an extra branch to the cenancestor.
 *
 * @author Diego Mallo
 * @version $Id: v 1.0 $
 */

public class CenancestorTreeLikelihood extends AbstractTreeLikelihood {
    private static final boolean DEBUG = false;
    private Parameter cenancestorHeight = null;
    private Parameter cenancestorBranch = null;
    private boolean branchRules=true;

    /**
     * Constructor.
     */
    public CenancestorTreeLikelihood(PatternList patternList,
                          TreeModel treeModel,
                          SiteModel siteModel,
                          CenancestorBranchRateModel branchRateModel,
                          TipStatesModel tipStatesModel,
                          Parameter cenancestorHeight,
                          Parameter cenancestorBranch,
                         // Parameter asStatistic,
                          boolean useAmbiguities,
                          boolean allowMissingTaxa,
                          boolean storePartials,
                          boolean forceJavaCore,
                          boolean forceRescaling) {

        super(CenancestorTreeLikelihoodParser.TREE_LIKELIHOOD, patternList, treeModel);

        this.storePartials = storePartials;
        nodeCount=treeModel.getNodeCount()+1;
        updateNode = new boolean[nodeCount];
        for (int i = 0; i < nodeCount; i++) {
            updateNode[i] = true;
        }

        try {
            
			final Logger logger = Logger.getLogger("dr.evomodel");
            
			this.siteModel = siteModel;
            addModel(siteModel);

            this.frequencyModel = siteModel.getFrequencyModel();
            addModel(frequencyModel);

            this.tipStatesModel = tipStatesModel;

            integrateAcrossCategories = siteModel.integrateAcrossCategories();

            this.categoryCount = siteModel.getCategoryCount();
            
            this.cenancestorHeight = cenancestorHeight;
            addVariable(cenancestorHeight);
            cenancestorHeight.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1)); 
            
            this.cenancestorBranch= cenancestorBranch;
            cenancestorBranch.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0,1));
            addVariable(cenancestorBranch);
            
			//Double-checking that the initial value of cenancestor-branch and otherwise making it appropriate (throwing a warning)
			double rootDepth=treeModel.getRootHeightParameter().getParameterValue(0);
			double initialCenancestorHeight=rootDepth+cenancestorBranch.getValue(0);
			double newCenancestorBranch=-1;
			
			if(initialCenancestorHeight<cenancestorHeight.getBounds().getLowerLimit(0)){
				newCenancestorBranch=cenancestorHeight.getBounds().getLowerLimit(0)-rootDepth;
			}
			if(initialCenancestorHeight>cenancestorHeight.getBounds().getUpperLimit(0) ){
				newCenancestorBranch=cenancestorHeight.getBounds().getUpperLimit(0)-rootDepth;
			}
			
			if(newCenancestorBranch!=-1){
				logger.info("The initial parameter value for the cenancestor branch length, "+cenancestorBranch.getValue(0)+" was not valid. Recalculated to: "+newCenancestorBranch);
				cenancestorBranch.setParameterValue(0, newCenancestorBranch);
			}

			branchRateModel.initCenancestor(cenancestorBranch);
            
            //if (asStatistic == cenancestorHeight){
            //	this.branchRules=true;
            //}
            
            //	if (branchRules==true){
            		updateCenancestorHeight(); //Trying to avoid improper initial values
            //	}
            // 	else{
            //		updateCenancestorBranch();
            //	}
        
            String coreName = "Java general";
            
            /**
             * TODO: Check if it is worth implementing other datatypes.
             */
            
            final DataType dataType = patternList.getDataType();
            
            if (dataType instanceof dr.evolution.datatype.TwoStates)
            {
            		coreName = "Java cenancestor binary";
            		cenancestorlikelihoodCore = new GeneralCenancestorLikelihoodCore(patternList.getStateCount());
            }
            else if (dataType instanceof dr.evolution.datatype.GeneralDataType) {
        			coreName = "Java cenancestor CNV";
        			cenancestorlikelihoodCore = new GeneralCenancestorLikelihoodCore(patternList.getStateCount());
            }

/*            if (integrateAcrossCategories) {

                final DataType dataType = patternList.getDataType();

                if (dataType instanceof dr.evolution.datatype.Nucleotides) {

                    if (!forceJavaCore && NativeNucleotideLikelihoodCore.isAvailable()) {
                        coreName = "native nucleotide";
                        likelihoodCore = new NativeNucleotideLikelihoodCore();
                    } else {
                        coreName = "Java nucleotide";
                        likelihoodCore = new NucleotideLikelihoodCore();
                    }

                } else if (dataType instanceof dr.evolution.datatype.AminoAcids) {
                    if (!forceJavaCore && NativeAminoAcidLikelihoodCore.isAvailable()) {
                        coreName = "native amino acid";
                        likelihoodCore = new NativeAminoAcidLikelihoodCore();
                    } else {
                        coreName = "Java amino acid";
                        likelihoodCore = new AminoAcidLikelihoodCore();
                    }

                    // The codon core was out of date and did nothing more than the general core...
                } else if (dataType instanceof dr.evolution.datatype.Codons) {
                    if (!forceJavaCore && NativeGeneralLikelihoodCore.isAvailable()) {
                        coreName = "native general";
                        likelihoodCore = new NativeGeneralLikelihoodCore(patternList.getStateCount());
                    } else {
                        coreName = "Java general";
                        likelihoodCore = new GeneralLikelihoodCore(patternList.getStateCount());
                    }
                    useAmbiguities = true;
                } else {
                    if (!forceJavaCore && NativeGeneralLikelihoodCore.isAvailable()) {
                        coreName = "native general";
                        likelihoodCore = new NativeGeneralLikelihoodCore(patternList.getStateCount());
                    } else {
                        	coreName = "Java general";
                        	likelihoodCore = new GeneralLikelihoodCore(patternList.getStateCount());
                    }
                }
            } else {
                likelihoodCore = new GeneralLikelihoodCore(patternList.getStateCount());
            }*/
            {
            		final String id = getId();
            		logger.info("TreeLikelihood(" + ((id != null) ? id : treeModel.getId()) + ") using " + coreName + " likelihood core");

            		logger.info("  " + (useAmbiguities ? "Using" : "Ignoring") + " ambiguities in tree likelihood.");
            		logger.info("  With " + patternList.getPatternCount() + " unique site patterns.");
            }

            if (branchRateModel != null) {
                this.branchRateModel = branchRateModel;
                logger.info("Branch rate model used: " + branchRateModel.getModelName());
            } else {
                this.branchRateModel = new DefaultCenancestorBranchRateModel();
            }
            addModel(this.branchRateModel);

            probabilities = new double[stateCount * stateCount];

            cenancestorlikelihoodCore.initialize(nodeCount, patternCount, categoryCount, integrateAcrossCategories);

            int extNodeCount = treeModel.getExternalNodeCount();
            int intNodeCount = treeModel.getInternalNodeCount();

            if (tipStatesModel != null) {
                tipStatesModel.setTree(treeModel);

                tipPartials = new double[patternCount * stateCount];

                for (int i = 0; i < extNodeCount; i++) {
                    // Find the id of tip i in the patternList
                    String id = treeModel.getTaxonId(i);
                    int index = patternList.getTaxonIndex(id);

                    if (index == -1) {
                        throw new TaxonList.MissingTaxonException("Taxon, " + id + ", in tree, " + treeModel.getId() +
                                ", is not found in patternList, " + patternList.getId());
                    }

                    tipStatesModel.setStates(patternList, index, i, id);
                    cenancestorlikelihoodCore.createNodePartials(i);
                }

                addModel(tipStatesModel);
            } else {
                for (int i = 0; i < extNodeCount; i++) {
                    // Find the id of tip i in the patternList
                    String id = treeModel.getTaxonId(i);
                    int index = patternList.getTaxonIndex(id);

                    if (index == -1) {
                        if (!allowMissingTaxa) {
                            throw new TaxonList.MissingTaxonException("Taxon, " + id + ", in tree, " + treeModel.getId() +
                                    ", is not found in patternList, " + patternList.getId());
                        }
                        if (useAmbiguities) {
                            setMissingPartials((LikelihoodCore)cenancestorlikelihoodCore, i);
                        } else {
                            setMissingStates((LikelihoodCore)cenancestorlikelihoodCore, i);
                        }
                    } else {
                        if (useAmbiguities) {
                            setPartials((LikelihoodCore)cenancestorlikelihoodCore, patternList, categoryCount, index, i);
                        } else {
                            setStates((LikelihoodCore)cenancestorlikelihoodCore, patternList, index, i);
                        }
                    }
                }
            }
            for (int i = 0; i <= intNodeCount; i++) { //Added one step for the cenancestor
                cenancestorlikelihoodCore.createNodePartials(extNodeCount + i);
            }
            
            
            if (forceRescaling) {
                cenancestorlikelihoodCore.setUseScaling(true);
                logger.info("  Forcing use of partials rescaling.");
            }

        } catch (TaxonList.MissingTaxonException mte) {
            throw new RuntimeException(mte.toString());
        }

        addStatistic(new SiteLikelihoodsStatistic());
    }

    public final CenancestorLikelihoodCore getLikelihoodCore() {
        return cenancestorlikelihoodCore;
    }
    
    //
    // Cenancestor-related functions
    //
    
    /**
     * set cenancestor date
     */
    public void setCenancestorHeight(double cen) {cenancestorHeight.setParameterValue(0, cen);}
    /**
     * @return cenancestor
     */
    public final double getCenancestorHeight() { return cenancestorHeight.getValue(0); }
    
    /**
     * set cenancestor_branch length
     */
    public void setCenancestorBranch(double cen) {cenancestorBranch.setParameterValue(0, cen);}
    /**
     * @return cenancestor
     */
    public final double getCenancestorBranch() { return cenancestorBranch.getValue(0); }
    
    private final void updateCenancestorBranch()
    {
    		cenancestorBranch.setParameterValueQuietly(0, getCenancestorHeight()-treeModel.getNodeHeight(treeModel.getRoot()));
    }
    
    private final void updateCenancestorHeight()
    {
    		cenancestorHeight.setParameterValueQuietly(0, treeModel.getNodeHeight(treeModel.getRoot())+getCenancestorBranch());
    }
    
    // **************************************************************
    // VariableListener IMPLEMENTATION
    // **************************************************************

    /**
     * This method is called whenever a parameter is changed.
     * <p/>
     * It is strongly recommended that the model component sets a "dirty" flag and does no
     * further calculations. Recalculation is typically done when the model component is asked for
     * some information that requires them. This mechanism is 'lazy' so that this method
     * can be safely called multiple times with minimal computational cost.
     */
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type)
    {
    		if (variable==this.cenancestorHeight)
    		{
    			updateCenancestorBranch();
    			updateNode(treeModel.getRoot());
    			fireModelChanged();
    		}
    		else if (variable==this.cenancestorBranch)
    		{
    			updateCenancestorHeight();
    			updateNode(treeModel.getRoot());
    			fireModelChanged();
    		}
    }
    
    
    // **************************************************************
    // ModelListener IMPLEMENTATION
    // **************************************************************

    /**
     * Handles model changed events from the submodels.
     */
    protected void handleModelChangedEvent(Model model, Object object, int index) {

        if (model == treeModel) {
            if (object instanceof TreeModel.TreeChangedEvent) {

            		if (((TreeModel.TreeChangedEvent) object).areAllInternalHeightsChanged()) {
                    updateAllNodes();
                    updateCenancestorHeight();
            		} else if (((TreeModel.TreeChangedEvent) object).isNodeChanged()) {
                    // If a node event occurs the node and its two child nodes
                    // are flagged for updating (this will result in everything
                    // above being updated as well. Node events occur when a node
                    // is added to a branch, removed from a branch or its height or
                    // rate changes.
                    updateNodeAndChildren(((TreeModel.TreeChangedEvent) object).getNode());
                    if(((TreeModel.TreeChangedEvent) object).getNode()==treeModel.getRoot()) {
                   		updateCenancestorHeight();
                    }

                } else if (((TreeModel.TreeChangedEvent) object).isTreeChanged()) {
                    // Full tree events result in a complete updating of the tree likelihood
                    updateAllNodes();
                    updateCenancestorHeight();
                } else {
                    // Other event types are ignored (probably trait changes).
                    //System.err.println("Another tree event has occured (possibly a trait change).");
                }
            }

        } else if (model == branchRateModel) {
            if (index == -1) {
                updateAllNodes();
            } else {
                if (DEBUG) {
                if (index >= treeModel.getNodeCount()) {
                    throw new IllegalArgumentException("Node index out of bounds");
                }
                }
                updateNode(treeModel.getNode(index));
            }

        } else if (model == frequencyModel) {

            updateAllNodes();

        } else if (model == tipStatesModel) {
        	if(object instanceof Taxon)
        	{
        		for(int i=0; i<treeModel.getNodeCount(); i++)
        			if(treeModel.getNodeTaxon(treeModel.getNode(i))!=null && treeModel.getNodeTaxon(treeModel.getNode(i)).getId().equalsIgnoreCase(((Taxon)object).getId()))
        				updateNode(treeModel.getNode(i));
        	}else
        		updateAllNodes();

        } else if (model instanceof SiteModel) {

            updateAllNodes();

        } else {

            throw new RuntimeException("Unknown componentChangedEvent");
        }

        super.handleModelChangedEvent(model, object, index);
    }

    // **************************************************************
    // Model IMPLEMENTATION
    // **************************************************************

    /**
     * Stores the additional state other than model components
     */
    protected void storeState() {

        if (storePartials) {
            cenancestorlikelihoodCore.storeState();
        }
        super.storeState();

    }

    /**
     * Restore the additional stored state
     */
    protected void restoreState() {

        if (storePartials) {
            cenancestorlikelihoodCore.restoreState();
        } else {
            updateAllNodes();
        }

        super.restoreState();

    }

    // **************************************************************
    // Likelihood IMPLEMENTATION
    // **************************************************************

    /**
     * Calculate the log likelihood of the current state.
     *
     * @return the log likelihood.
     */
    protected double calculateLogLikelihood() {

        if (patternLogLikelihoods == null) {
            patternLogLikelihoods = new double[patternCount];
        }

        if (!integrateAcrossCategories) {
            if (siteCategories == null) {
                siteCategories = new int[patternCount];
            }
            for (int i = 0; i < patternCount; i++) {
                siteCategories[i] = siteModel.getCategoryOfSite(i);
            }
        }

        if (tipStatesModel != null) {
            int extNodeCount = treeModel.getExternalNodeCount();
            for (int index = 0; index < extNodeCount; index++) {
                if (updateNode[index]) {
                    cenancestorlikelihoodCore.setNodePartialsForUpdate(index);
                    tipStatesModel.getTipPartials(index, tipPartials);
                    cenancestorlikelihoodCore.setCurrentNodePartials(index, tipPartials);
                }
            }
        }


        final NodeRef root = treeModel.getRoot();
        traverse(treeModel, root);

        double logL = 0.0;
        double ascertainmentCorrection = getAscertainmentCorrection(patternLogLikelihoods);
        for (int i = 0; i < patternCount; i++) {
            logL += (patternLogLikelihoods[i] - ascertainmentCorrection) * patternWeights[i];
        }

        if (logL == Double.NEGATIVE_INFINITY) {
            Logger.getLogger("dr.evomodel").info("TreeLikelihood, " + this.getId() + ", turning on partial likelihood scaling to avoid precision loss");

            // We probably had an underflow... turn on scaling
            cenancestorlikelihoodCore.setUseScaling(true);

            // and try again...
            updateAllNodes();
            updateAllPatterns();
            traverse(treeModel, root);

            logL = 0.0;
            ascertainmentCorrection = getAscertainmentCorrection(patternLogLikelihoods);
            for (int i = 0; i < patternCount; i++) {
                logL += (patternLogLikelihoods[i] - ascertainmentCorrection) * patternWeights[i];
            }
        }

        //********************************************************************
        // after traverse all nodes and patterns have been updated --
        //so change flags to reflect this.
        for (int i = 0; i < nodeCount; i++) {
            updateNode[i] = false;
        }
        //********************************************************************

        return logL;
    }

    public double[] getPatternLogLikelihoods() {
        getLogLikelihood(); // Ensure likelihood is up-to-date
        double ascertainmentCorrection = getAscertainmentCorrection(patternLogLikelihoods);
        double[] out = new double[patternCount];
        for (int i = 0; i < patternCount; i++) {
            if (patternWeights[i] > 0) {
                out[i] = (patternLogLikelihoods[i] - ascertainmentCorrection) * patternWeights[i];
            } else {
                out[i] = Double.NEGATIVE_INFINITY;
            }
        }
        return out;
    }

    /* Calculate ascertainment correction if working off of AscertainedSitePatterns
    @param patternLogProbs log pattern probabilities
    @return the log total probability for a pattern.
    */
    protected double getAscertainmentCorrection(double[] patternLogProbs) {
        if (patternList instanceof AscertainedSitePatterns) {
            return ((AscertainedSitePatterns) patternList).getAscertainmentCorrection(patternLogProbs);
        } else {
            return 0.0;
        }
    }

    /**
     * Check whether the scaling is still required. If the sum of all the logScalingFactors
     * is zero then we simply turn off the useScaling flag. This will speed up the likelihood
     * calculations when scaling is not required.
     */
    public void checkScaling() {
//	    if (useScaling) {
//	        if (scalingCheckCount % 1000 == 0) {
//	            double totalScalingFactor = 0.0;
//	            for (int i = 0; i < nodeCount; i++) {
//	                for (int j = 0; j < patternCount; j++) {
//	                    totalScalingFactor += scalingFactors[currentPartialsIndices[i]][i][j];
//	                }
//	            }
//	            useScaling = totalScalingFactor < 0.0;
//	            Logger.getLogger("dr.evomodel").info("LikelihoodCore total log scaling factor: " + totalScalingFactor);
//	            if (!useScaling) {
//	                Logger.getLogger("dr.evomodel").info("LikelihoodCore scaling turned off.");
//	            }
//	        }
//	        scalingCheckCount++;
//	    }
    }


    /**
     * Traverse the tree calculating partial likelihoods.
     *
     * @return whether the partials for this node were recalculated.
     */
    protected boolean traverse(Tree tree, NodeRef node) {

        boolean update = false;
        boolean rootUpdated = false;

        int nodeNum = node.getNumber();

        NodeRef parent = tree.getParent(node);

        // First update the transition probability matrix(ices) for this branch if it is a normal branch
        if (parent != null && updateNode[nodeNum]) {

            final double branchRate = branchRateModel.getBranchRate(tree, node);

            // Get the operational time of the branch
            final double branchTime = branchRate * (tree.getNodeHeight(parent) - tree.getNodeHeight(node));

            if (branchTime < 0.0) {
                throw new RuntimeException("Negative branch length: " + branchTime);
            }

            cenancestorlikelihoodCore.setNodeMatrixForUpdate(nodeNum);

            for (int i = 0; i < categoryCount; i++) {

                double branchLength = siteModel.getRateForCategory(i) * branchTime;
                siteModel.getSubstitutionModel().getTransitionProbabilities(branchLength, probabilities);
                cenancestorlikelihoodCore.setNodeMatrix(nodeNum, i, probabilities);
            }

            update = true;
        }
        else if (parent == null && cenancestorHeight != null && updateNode[nodeNum]) //The root has to be updated
        {
    			// First update the transition probability matrix(ices) for the root-cenancestor fake branch
    			rootUpdated=true;
    			// Get the operational time of the fake branch from the root to the cenancestor 
    			double rootHeight = treeModel.getNodeHeight(treeModel.getRoot());
    			double branchRate = branchRateModel.getBranchRate(rootHeight, getCenancestorHeight()); //TODO: Could this be easily improved? I would have to adapt the tree structure and abstact tree likelihood
    			double branchTime = branchRate * getCenancestorBranch() ; //TODO: Could this be easily improved? The same as before

    		
    			for (int i = 0; i < categoryCount; i++) 
    			{
    				double branchLength = siteModel.getRateForCategory(i) * branchTime;
    				siteModel.getSubstitutionModel().getTransitionProbabilities(branchLength, probabilities);
    				cenancestorlikelihoodCore.setNodeMatrix(nodeNum, i, probabilities);

    			}
        }

        // If the node is internal, update the partial likelihoods.
        if (!tree.isExternal(node)) {

            // Traverse down the two child nodes
            NodeRef child1 = tree.getChild(node, 0);
            final boolean update1 = traverse(tree, child1);

            NodeRef child2 = tree.getChild(node, 1);
            final boolean update2 = traverse(tree, child2);

            // If either child node was updated then update this node too
            if (update1 || update2 || rootUpdated) {

                if (update1 || update2) {
            			final int childNum1 = child1.getNumber();
            			final int childNum2 = child2.getNumber();

            			cenancestorlikelihoodCore.setNodePartialsForUpdate(nodeNum);

            			if (integrateAcrossCategories) {
            				cenancestorlikelihoodCore.calculatePartials(childNum1, childNum2, nodeNum);
            			} else {
            				cenancestorlikelihoodCore.calculatePartials(childNum1, childNum2, nodeNum, siteCategories);
            			}

            			if (COUNT_TOTAL_OPERATIONS) {
            				totalOperationCount ++;
            			}
                }
                
                if (parent == null) {
                    // No parent this is the root of the tree
                	
                    double[] partials;
                    int nodeNumCenan = getCenancestorIndex();
                    
                    if(cenancestorHeight != null) {
                    		if (update1 || update2 || rootUpdated) {
                        		// Calculate the partials at the cenancestor. The transition matrix of the root was calculated before.
                        		cenancestorlikelihoodCore.setNodePartialsForUpdate(nodeNumCenan);

                        		if (integrateAcrossCategories) {
                        			cenancestorlikelihoodCore.calculatePartials(nodeNum, nodeNumCenan);
                        		} else {
                        			cenancestorlikelihoodCore.calculatePartials(nodeNum, nodeNumCenan, siteCategories);
                        		}
                    		}
                    		
                    		partials = getCenancestorPartials();
                    }
                    else { //Using the cenancestor model without cenancestor date. It assumes that the root of the tree is the cenancestor. Not tested. It shouldn't be normally used either.
                    		partials = getRootPartials();
                    }
                    
                    	// calculate the pattern likelihoods
                    double[] frequencies = frequencyModel.getFrequencies();
                    cenancestorlikelihoodCore.calculateLogLikelihoods(partials, frequencies, patternLogLikelihoods);
                }

                update = true;
            }
        }

        return update;

    }

    public final double[] getRootPartials() {
        if (rootPartials == null) {
            rootPartials = new double[patternCount * stateCount];
        }

        int nodeNum = treeModel.getRoot().getNumber();
        if (integrateAcrossCategories) {

            // moved this call to here, because non-integrating siteModels don't need to support it - AD
            double[] proportions = siteModel.getCategoryProportions();
            		cenancestorlikelihoodCore.integratePartials(nodeNum, proportions, rootPartials);
        } else {
        			cenancestorlikelihoodCore.getPartials(nodeNum, rootPartials);
        }

        return rootPartials;
    }
    
    public final double[] getCenancestorPartials() {
        if (cenancestorPartials == null) {
            cenancestorPartials = new double[patternCount * stateCount];
        }

        int nodeNum = getCenancestorIndex();
        if (integrateAcrossCategories) {

            // moved this call to here, because non-integrating siteModels don't need to support it - AD
            double[] proportions = siteModel.getCategoryProportions();
            		cenancestorlikelihoodCore.integratePartials(nodeNum, proportions, cenancestorPartials);
        } else {
        			cenancestorlikelihoodCore.getPartials(nodeNum, cenancestorPartials);
        }

        return cenancestorPartials;
    }
 
    
    /**
     * Return the index of the fake cenancestor node to modify its associated matrices and partials.
     * @return Index of fake cenancestor node
     */
    public int getCenancestorIndex() {
    		return nodeCount-1;
    }

    /**
     * the root partial likelihoods (a temporary array that is used
     * to fetch the partials - it should not be examined directly -
     * use getRootPartials() instead).
     */
    private double[] rootPartials = null;
    private double[] cenancestorPartials = null;

    public class SiteLikelihoodsStatistic extends Statistic.Abstract {

        public SiteLikelihoodsStatistic() {
            super("siteLikelihoods");
        }

        public int getDimension() {
            if (patternList instanceof SitePatterns) {
                return ((SitePatterns)patternList).getSiteCount();
            } else {
                return patternList.getPatternCount();
            }
        }

        public String getDimensionName(int dim) {
            return getTreeModel().getId() + "site-" + dim;
        }

        public double getStatisticValue(int i) {

            if (patternList instanceof SitePatterns) {
                int index = ((SitePatterns)patternList).getPatternIndex(i);
                if( index >= 0 ) {
                    return patternLogLikelihoods[index] / patternWeights[index];
                } else {
                    return 0.0;
                }
            } else {
                return patternList.getPatternCount();
            }
        }
    }

    // **************************************************************
    // INSTANCE VARIABLES
    // **************************************************************

    /**
     * the frequency model for these sites
     */
    protected final FrequencyModel frequencyModel;

    /**
     * the site model for these sites
     */
    protected final SiteModel siteModel;

    /**
     * the branch rate model
     */
    protected final CenancestorBranchRateModel branchRateModel;

    /**
     * the tip partials model
     */
    private final TipStatesModel tipStatesModel;

    private final boolean storePartials;

    protected final boolean integrateAcrossCategories;

    /**
     * the categories for each site
     */
    protected int[] siteCategories = null;


    /**
     * the pattern likelihoods
     */
    protected double[] patternLogLikelihoods = null;

    /**
     * the number of rate categories
     */
    protected int categoryCount;

    /**
     * an array used to transfer transition probabilities
     */
    protected double[] probabilities;


    /**
     * an array used to transfer tip partials
     */
    protected double[] tipPartials;

    /**
     * the LikelihoodCore
     */
    protected CenancestorLikelihoodCore cenancestorlikelihoodCore;
}
