package ga_complex_planning;

import java.util.HashSet;
import java.util.Random;
import java.util.Set;

/**
 * @ClassName ConstraintChecker
 * @Description 检查基因组是否超出约束
 * @Author fwt
 * @Date 2024/1/25 4:35 下午
 * @Version 1.0
 **/
/*
public class ConstraintChecker {
    private static int i = 0;

    public static boolean isGoodChromosome(Chromosome chromosome) {
        // 频段约束
        if (!checkFrequencyRangeConstraint(chromosome.getFrequencies(), chromosome.getF_min(), chromosome.getF_max())) {
            System.out.println("频率出界");
            return false;
        }

        // 频点范围约束
        if (!checkFrequencyPointRangeConstraint(chromosome.getFrequencies(), chromosome.getF_disable_min(), chromosome.getF_disable_max())) {
            System.out.println("频率禁忌");
            return false;

        }

        // 带宽上限约束
        if (!checkBandwidthLimitConstraint(chromosome.getBandwidths(), chromosome.getB_max())) {
//            System.out.println("带宽出界");
            return false;
        }

        return true;
    }

    private static boolean checkFrequencyRangeConstraint(int[] frequencies, double fMin, double fMax) {
        for (double frequency : frequencies) {
            if (frequency < fMin || frequency > fMax) {

                return false;
            }
        }
        return true;
    }

    private static boolean checkFrequencyPointRangeConstraint(int[] frequencies, double fRangeMin, double fRangeMax) {
        for (double frequency : frequencies) {
            if (frequency >= fRangeMin && frequency <= fRangeMax) {

                return false; // 如果有任何一个频率在范围内，返回 false，表示不符合约束
            }
        }
        return true; // 所有频率都不在范围内，返回 true，表示满足约束
    }

    private static boolean checkBandwidthLimitConstraint(int[] bandwidths, double bMax) {
        i++;
        double totalBandwidth = 0.0;
        for (double bandwidth : bandwidths) {
//            System.out.println(bandwidth);
            totalBandwidth += bandwidth;
        }
//        System.out.println("第"+i+"个基因总带宽为"+totalBandwidth);
        return totalBandwidth <= bMax;
    }

}
*/


public class ConstraintChecker {

    public static boolean isGoodChromosome(Chromosome chromosome) {
        boolean isGood = true;

        // 频段约束
        if (!checkFrequencyRangeConstraint(chromosome.getFrequencies(), chromosome.getF_min(), chromosome.getF_max())) {
//            System.out.println("频率出界");
            isGood = false;
        }

        // 频点范围约束
        if (!checkFrequencyPointRangeConstraint(chromosome.getFrequencies(), chromosome.getF_disable_min(), chromosome.getF_disable_max())) {
//            System.out.println("频率禁忌");
            isGood = false;
        }

        // 频率唯一性约束
        if (!checkFrequencyUniqueness(chromosome.getFrequencies())) {
//            System.out.println("频率不唯一");
            isGood = false;
        }

        // 带宽上限约束
        if (!checkBandwidthLimitConstraint(chromosome.getBandwidths(), chromosome.getB_max())) {
//            System.out.println("带宽出界");
            isGood = false;
        }

        return isGood;
    }

    public static void enforceConstraints(Chromosome chromosome) {
        // 修正频率范围
        fixFrequencyRange(chromosome.getFrequencies(), chromosome.getF_min(), chromosome.getF_max());

        // 修正禁用频率范围
        fixFrequencyPointRange(chromosome.getFrequencies(), chromosome.getF_disable_min(), chromosome.getF_disable_max());

        // 修正频率唯一性
        ensureFrequencyUniqueness(chromosome.getFrequencies(), chromosome.getF_min(), chromosome.getF_max());

        // 修正带宽限制
        fixBandwidthLimit(chromosome.getBandwidths(), chromosome.getB_max());
    }

    // 修正频率范围
    private static void fixFrequencyRange(int[] frequencies, double fMin, double fMax) {
        for (int i = 0; i < frequencies.length; i++) {
            if (frequencies[i] < fMin || frequencies[i] > fMax) {
                frequencies[i] = (int) Math.max(fMin, Math.min(fMax, frequencies[i]));
            }
        }
    }

    // 修正禁用频率范围
    private static void fixFrequencyPointRange(int[] frequencies, double fRangeMin, double fRangeMax) {
        for (int i = 0; i < frequencies.length; i++) {
            if (frequencies[i] >= fRangeMin && frequencies[i] <= fRangeMax) {
                frequencies[i] = (int) fRangeMin - 1;  // 将频率调整到禁用范围外
            }
        }
    }

    private static void ensureFrequencyUniqueness(int[] frequencies, double fMin, double fMax) {
        Set<Integer> frequencySet = new HashSet<>();
        Random random = new Random();

        for (int i = 0; i < frequencies.length; i++) {
            while (!frequencySet.add(frequencies[i])) {
                // 如果频率不唯一，重新生成一个频率
                frequencies[i] = (int) (fMin + random.nextDouble() * (fMax - fMin));
            }
        }
    }


    // 检查频率的唯一性
    private static boolean checkFrequencyUniqueness(int[] frequencies) {
        Set<Integer> frequencySet = new HashSet<>();
        for (int frequency : frequencies) {
            if (!frequencySet.add(frequency)) {
                return false;  // 如果添加失败，说明频率重复
            }
        }
        return true;
    }

    // 修正带宽限制
    private static void fixBandwidthLimit(int[] bandwidths, double bMax) {
        double totalBandwidth = 0.0;
        for (int bandwidth : bandwidths) {
            totalBandwidth += bandwidth;
        }

        if (totalBandwidth > bMax) {
            double scale = bMax / totalBandwidth;
            for (int i = 0; i < bandwidths.length; i++) {
                bandwidths[i] *= scale;
            }
        }
    }

    private static boolean checkFrequencyRangeConstraint(int[] frequencies, double fMin, double fMax) {
        for (double frequency : frequencies) {
            if (frequency < fMin || frequency > fMax) {

                return false;
            }
        }
        return true;
    }

    private static boolean checkFrequencyPointRangeConstraint(int[] frequencies, double fRangeMin, double fRangeMax) {
        for (double frequency : frequencies) {
            if (frequency >= fRangeMin && frequency <= fRangeMax) {

                return false; // 如果有任何一个频率在范围内，返回 false，表示不符合约束
            }
        }
        return true; // 所有频率都不在范围内，返回 true，表示满足约束
    }

    private static boolean checkBandwidthLimitConstraint(int[] bandwidths, double bMax) {
//        i++;
        double totalBandwidth = 0.0;
        for (double bandwidth : bandwidths) {
//            System.out.println(bandwidth);
            totalBandwidth += bandwidth;
        }
//        System.out.println("第"+i+"个基因总带宽为"+totalBandwidth);
        return totalBandwidth <= bMax;
    }

    public boolean areChannelsMutuallyExclusive(int[] frequencies, int[] bandwidths) {
        for (int i = 0; i < frequencies.length; i++) {
            for (int j = i + 1; j < frequencies.length; j++) {
                double channelStart1 = frequencies[i] - bandwidths[i] / 2.0;
                double channelEnd1 = frequencies[i] + bandwidths[i] / 2.0;
                double channelStart2 = frequencies[j] - bandwidths[j] / 2.0;
                double channelEnd2 = frequencies[j] + bandwidths[j] / 2.0;

                int tolerance = 1500;
                if (!(channelEnd1 + tolerance <= channelStart2 || channelEnd2 + tolerance <= channelStart1)) {
                    return false;
                }
            }
        }
        return true;
    }
}

