/*
 * Copyright (C) 2020. Tim Vaughan
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

package pitchfork.models;

import beast.base.core.Input;
import beast.base.evolution.tree.IntervalType;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeDistribution;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import beast.base.util.Binomial;

import java.math.BigDecimal;

public class BetaSkylineDistributionVersion3 extends AbstractBetaSkylineDistribution {

    public BetaSkylineDistributionVersion3() {
        super();
    }


    // This version updates N after a certain time
    @Override
    public double calculateLogP() {
        logP = 0.0;
        double T = tree.getRoot().getHeight();
        //System.out.print("New line");
        //System.out.print(T);

        double t = 0;
        double N = skylinePopulationsInput.get().getValue(0);
        //int seenCoalescentEvents = 0;
        double dt_track = 0;                                        // this variable is first dt after every round and then t_update is deducted until dt_track < t_update
        double t_update = T/populationSizes.getDimension();                                     //this variable defines the length of the time intervals
        int counter = 0;                                             // counter for groups
        double t_rest = 0;      //this variable is needed that we don't use dt_track after the while loop AND again in the next for loop round

        for (int i = 0; i < collapsedTreeIntervals.getIntervalCount(); i++) {
            // Get interval details
            double dt = collapsedTreeIntervals.getInterval(i);
            dt_track += dt;

            int n = collapsedTreeIntervals.getLineageCount(i);
            t += dt;
            //System.out.print(t + "\t");

            while(dt_track > t_update && counter < skylinePopulationsInput.get().getDimension()-1 ){
                logP += Math.log(betaCoalescentModel.getTotalCoalRate(n))-betaCoalescentModel.getTotalCoalRate(n)*(t_update-t_rest)/N;  // here we calculate the logP of the waiting time for the pre-defined interval t_update as long dt_track is bigger than t_update
                dt_track -= t_update;   // in the first round of the while loop we need to update ealier because there is still some dt_track to deduct from the last round
                t_rest = 0;           //now t_rest should be set to 0 again
                //if(dt_track>=2*t_update | t<T) {
                counter += 1;
                N = skylinePopulationsInput.get().getValue(counter);
                //}
            }

            t_rest = dt_track;
            logP += -betaCoalescentModel.getTotalCoalRate(n)*t_rest/N;      // here we calculate the logP for the time interval which is smaller or equal than t_update

            if (collapsedTreeIntervals.getIntervalType(i) == IntervalType.COALESCENT) {
                // Beta-coalescent event contribution

                int k = collapsedTreeIntervals.getCoalescentEvents(i)+1;
                logP += betaCoalescentModel.getLogLambda(n, k) + Binomial.logChoose(n, k) - Math.log(N) - Math.log(betaCoalescentModel.getTotalCoalRate(n));
            }
        }
        return logP;
    }

    @Override
    protected boolean requiresRecalculation() {
        return true;
    }

    @Override
    double[] getPopSizes(int gridSize) {
        double[] popSizes = new double[gridSize];
        double T = tree.getRoot().getHeight();

        int group = 0;
        double N = skylinePopulationsInput.get().getValue(0);
        int i = 0;
        double t_total = 0;
        double t_update = T/populationSizes.getDimension();
        double dt_track = 0;

        for (int gridIdx = 0; gridIdx < gridSize; gridIdx++) {
            double t = gridIdx * T / (gridSize - 1);

            // TODO: update N as needed
            while((t_total+collapsedTreeIntervals.getInterval(i))<t){

                double dt = collapsedTreeIntervals.getInterval(i);

                // Increment time
                t_total += dt;
                dt_track += dt;
                while(dt_track > t_update && group < skylinePopulationsInput.get().getDimension()-1){
                    dt_track -= t_update;
                    //if(dt_track>=2*t_update | t_total<T) {
                    group += 1;
                    N = skylinePopulationsInput.get().getValue(group);
                    //}
                }

                if (i == (collapsedTreeIntervals.getIntervalCount()-1)) {  //this if statement ends while loop when we reach the root of the tree
                    break;
                }
                i++;

            }

            popSizes[gridIdx] = N;
        }

        return popSizes;
    }
}