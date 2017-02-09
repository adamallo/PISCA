/*
 * RandomLocalClockModelCenancestor.java
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

package pluginSCA;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.tree.TreeParameterModel;
import dr.evomodel.tree.randomlocalmodel.RandomLocalTreeVariable;
import dr.evomodelxml.branchratemodel.RandomLocalClockModelParser;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Statistic;
import dr.inference.model.Variable;

import java.util.logging.Logger;

/**
 * A model of rates evolving along a tree, such that at each node the rate may change or not depending on an
 * indicator which chooses whether the parent rate is inherited or a new rate begins on the branch above the node.
 * This model is implemented using stochastic variable selection. Adapted for the cenancestor tree model.
 *
 * @author Alexei Drummond
 * @author Andrew Rambaut
 * @author Diego Mallo
 * 
 */
public class RandomLocalClockModelCenancestor extends AbstractCenancestorBranchRateModel
        implements RandomLocalTreeVariable {

    public RandomLocalClockModelCenancestor(TreeModel treeModel,
                                 Parameter meanRateParameter,
                                 Parameter rateIndicatorParameter,
                                 Parameter ratesParameter,
                                 boolean ratesAreMultipliers) {

        super(RandomLocalClockModelParser.LOCAL_BRANCH_RATES);

        this.ratesAreMultipliers = ratesAreMultipliers;

        indicators = new TreeParameterModel(treeModel, rateIndicatorParameter, false);
        rates = new TreeParameterModel(treeModel, ratesParameter, false);

        rateIndicatorParameter.addBounds(new Parameter.DefaultBounds(1, 0, rateIndicatorParameter.getDimension()));
        ratesParameter.addBounds(new Parameter.DefaultBounds(Double.MAX_VALUE, 0, ratesParameter.getDimension()));

        for (int i = 0; i < rateIndicatorParameter.getDimension(); i++) {
            rateIndicatorParameter.setParameterValue(i, 0.0);
            ratesParameter.setParameterValue(i, 1.0);
        }

        this.meanRateParameter = meanRateParameter;

        addModel(treeModel);
        this.treeModel = treeModel;

        addModel(indicators);
        addModel(rates);
        if (meanRateParameter != null) addVariable(meanRateParameter);

        unscaledBranchRates = new double[treeModel.getNodeCount()];

        Logger.getLogger("dr.evomodel").info("  indicator parameter name is '" + rateIndicatorParameter.getId() + "'");

        initialScaleFactor();
    }

    /**
     * @param tree the tree
     * @param node the node to retrieve the variable of
     * @return the raw real-valued variable at this node
     */
    public final double getVariable(Tree tree, NodeRef node) {
        return rates.getNodeValue(tree, node);
    }

    /**
     * @param tree the tree
     * @param node the node
     * @return true of the variable at this node is included in function, thus representing a change in the
     *         function looking down the tree.
     */
    public final boolean isVariableSelected(Tree tree, NodeRef node) {
        return indicators.getNodeValue(tree, node) > 0.5;
    }

    public void handleModelChangedEvent(Model model, Object object, int index) {
        recalculationNeeded = true;
        fireModelChanged();
    }

    protected final void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        recalculationNeeded = true;
        fireModelChanged();
    }

    protected void storeState() {
    }

    protected void restoreState() {
        recalculateScaleFactor();
    }

    protected void acceptState() {
    }

    public double getBranchRate(final Tree tree, final NodeRef node) {
        if (recalculationNeeded) {
            recalculateScaleFactor();
            recalculationNeeded = false;
        }
        return unscaledBranchRates[node.getNumber()] * scaleFactor;
    }
    
    //The scaled branch rate for the MRCA is 1*scaleFactor
    public double getBranchRate(double mrcaHeight, double cenHeight) {
        if (recalculationNeeded) {
            recalculateScaleFactor();
            recalculationNeeded = false;
        }
        return scaleFactor;
    }

    private void calculateUnscaledBranchRates(TreeModel tree) {
        cubr(tree, tree.getRoot(), 1.0);//Unscaled branch rate of the root is 1 by definition
    }

    /**
     * This is a recursive function that does the work of
     * calculating the unscaled branch rates across the tree
     * taking into account the indicator variables.
     *
     * @param tree the tree
     * @param node the node
     * @param rate the rate of the parent node
     */
    private void cubr(TreeModel tree, NodeRef node, double rate) {

        int nodeNumber = node.getNumber();

        if (!tree.isRoot(node)) {
            if (isVariableSelected(tree, node)) {
                if (ratesAreMultipliers) {
                    rate *= getVariable(tree, node);
                } else {
                    rate = getVariable(tree, node);
                }
            }
        }
        unscaledBranchRates[nodeNumber] = rate;

        int childCount = tree.getChildCount(node);
        for (int i = 0; i < childCount; i++) {
            cubr(tree, tree.getChild(node, i), rate);
        }
    }

    private void recalculateScaleFactor() {

        calculateUnscaledBranchRates(treeModel);

        double timeTotal = 0.0;
        double branchTotal = 0.0;

        for (int i = 0; i < treeModel.getNodeCount(); i++) {
            NodeRef node = treeModel.getNode(i);
            if (!treeModel.isRoot(node)) {

                double branchInTime =
                        treeModel.getNodeHeight(treeModel.getParent(node)) -
                                treeModel.getNodeHeight(node);

                double branchLength = branchInTime * unscaledBranchRates[node.getNumber()];

                timeTotal += branchInTime;
                branchTotal += branchLength;
            }
            else { //For MRCA-LUCA branch
                double branchInTime =cenancestorBranch.getStatisticValue(0);

                double branchLength = branchInTime * unscaledBranchRates[node.getNumber()]; //This branch has rate0

                timeTotal += branchInTime;
                branchTotal += branchLength;
            }
        }

        scaleFactor = timeTotal / branchTotal;

        if (meanRateParameter != null)
            scaleFactor *= meanRateParameter.getParameterValue(0);
    }
    
    private void initialScaleFactor() { //The initial calculation of the scale factor lacks the cenancestor (and therefore it is incorrect, but it is just for initialization).

        calculateUnscaledBranchRates(treeModel);

        double timeTotal = 0.0;
        double branchTotal = 0.0;

        for (int i = 0; i < treeModel.getNodeCount(); i++) {
            NodeRef node = treeModel.getNode(i);
            if (!treeModel.isRoot(node)) {

                double branchInTime =
                        treeModel.getNodeHeight(treeModel.getParent(node)) -
                                treeModel.getNodeHeight(node);

                double branchLength = branchInTime * unscaledBranchRates[node.getNumber()];

                timeTotal += branchInTime;
                branchTotal += branchLength;
            }
            
        }

        scaleFactor = timeTotal / branchTotal;

        if (meanRateParameter != null)
            scaleFactor *= meanRateParameter.getParameterValue(0);
        
        recalculationNeeded=true;
    }
    
    public void initCenancestor(Parameter cenancestorBranch)
    {
    		this.cenancestorBranch=cenancestorBranch;
    		addVariable(cenancestorBranch); //This would actually be part of the tree model. If it changes, we need to recalculate the rates, and then the likelihood.
    }

    public Parameter getCenancestor() {
		return this.cenancestorBranch;
    }

    //Cenancestor-related variables
    private Parameter cenancestorBranch = null;
    
    // the scale factor necessary to maintain the mean rate
    private double scaleFactor;

    // the tree model
    private TreeModel treeModel;

    // true if the rate variables are treated as relative
    // to the parent rate rather than absolute rates
    private boolean ratesAreMultipliers = false;

    // the unscaled rates of each branch, taking into account the indicators
    private double[] unscaledBranchRates;

    // the mean rate across all the tree, if null then mean rate is scaled to 1.0
    private Parameter meanRateParameter;

    private TreeParameterModel indicators;
    private TreeParameterModel rates;

    boolean recalculationNeeded = true;
}
