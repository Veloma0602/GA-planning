package ga_complex_planning;

/**
 * @author fwt
 * @date 2024/8/17
 * @Description 计算综合的交叉频道干扰惩罚项
 */
public class CombinedPenaltyCalculator {

    private double penaltyMin;
    private double penaltyMax;

    public CombinedPenaltyCalculator() {
        penaltyMin = Double.MAX_VALUE;
        penaltyMax = Double.MIN_VALUE;
    }

    // 在初步探索时更新 P_min 和 P_max
    public void updateMinMax(int[] frequencies, int[] bandwidths) {
        double penalty = calculateRawPenalty(frequencies, bandwidths);
        penaltyMin = Math.min(penaltyMin, penalty);
        penaltyMax = Math.max(penaltyMax, penalty);
    }

    // 计算未归一化的综合惩罚项
    public double calculateRawPenalty(int[] frequencies, int[] bandwidths) {
        double penalty = 0.0;
        double epsilon = 1e-6;

        for (int i = 0; i < frequencies.length; i++) {
            for (int j = i + 1; j < frequencies.length; j++) {
                double freqDifference = Math.abs(frequencies[i] - frequencies[j]);
                double channelStart1 = frequencies[i] - bandwidths[i] / 2.0;
                double channelEnd1 = frequencies[i] + bandwidths[i] / 2.0;
                double channelStart2 = frequencies[j] - bandwidths[j] / 2.0;
                double channelEnd2 = frequencies[j] + bandwidths[j] / 2.0;

                double overlap = Math.max(0, Math.min(channelEnd1, channelEnd2) - Math.max(channelStart1, channelStart2));
                penalty += (1.0 / (freqDifference + epsilon)) * overlap;
            }
        }
        return penalty;
    }

    // 计算归一化的综合惩罚项
    public double calculateCombinedPenalty(int[] frequencies, int[] bandwidths, double P_min, double P_max) {
        double rawPenalty = calculateRawPenalty(frequencies, bandwidths);
        return (rawPenalty - P_min) / (P_max - P_min + 1e-6);  // 防止分母为零
    }
}

