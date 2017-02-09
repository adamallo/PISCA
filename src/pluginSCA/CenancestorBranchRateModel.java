/*
 * CenancestorBranchRateModel.java
 *
 * Modified by DM from BranchRateModel.java
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

import dr.evomodel.branchratemodel.BranchRateModel;
import dr.inference.model.Parameter;

/**
 * Date: Dec 13, 2004
 * Time: 1:59:24 PM
 *
 * @author Alexei Drummond
 * @version $Id: BranchRateModel.java,v 1.4 2005/05/24 20:25:57 rambaut Exp $
 */
public interface CenancestorBranchRateModel extends BranchRateModel {
    public static final String BRANCH_RATES = "branchRates";
    public static final String RATE = "rate";

    // This is inherited from BranchRates:
    //double getBranchRate(Tree tree, NodeRef node);
    double getBranchRate(double mrcaHeight, double cenHeight);
    
    //This connects the cenancestor parameter with BranchRateModels, since they cannot get it from the tree
    void initCenancestor(Parameter cenancestorBranch);
    
    //To spread the cenancestor to other classes using CenancestorBranchRateModel: Important, before all models are conected to the likelihood model
    //this my return null. All objects calling this should have a flag and keep trying to get the cenancestor until it is not null.
    Parameter getCenancestor();
}
