/*
 * CenancestorLikelihoodCore.java
 *
 * By Diego Mallo
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
import dr.evomodel.treelikelihood.LikelihoodCore;

/**
 * CenancestorLikelihoodCore - An interface extending the core likelihood functions for trees with cenancestor
 *
 * @author Diego Mallo
 */

public interface CenancestorLikelihoodCore extends LikelihoodCore {
    
    /**
     * Calculates partial likelihoods at a node with only one child
     *
     * @param nodeIndex1 the 'child 1' node
     * @param nodeIndex3 the 'parent' node
     */
    void calculatePartials(int nodeIndex1, int nodeIndex3);

    /**
     * Calculates partial likelihoods at a node with only one child using a matrixMap.
     *
     * @param nodeIndex1 the 'child 1' node
     * @param nodeIndex3 the 'parent' node
     * @param matrixMap  a map of which matrix to use for each pattern
     */
    void calculatePartials(int nodeIndex1,int nodeIndex3, int[] matrixMap);

}
