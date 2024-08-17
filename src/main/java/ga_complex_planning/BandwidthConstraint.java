package ga_complex_planning;

import util.PropertyUtil;

import java.util.Properties;

/**
 * @ClassName BandwidthConstraint
 * @Description 带宽限制函数
 * @Author fwt
 * @Date 2024/1/24 9:35 下午
 * @Version 1.0
 **/

public class BandwidthConstraint {

    // 定值
    private double lambda;
    private double b_max;
    private double F_bandwidth_max;  // 新增：带宽惩罚的最大值

    private Properties planningInfoPro = PropertyUtil.getProperty("planning_info");

    private int numOfPlatforms = Integer.valueOf(planningInfoPro.getProperty("numOfPlatforms"));


    // 构造函数
    public BandwidthConstraint(Chromosome chromosome) {
        this.lambda = chromosome.getLambda();
        this.b_max = chromosome.getB_max();

        this.F_bandwidth_max = lambda * numOfPlatforms * Math.pow(b_max, 2);
    }

    // 计算带宽限制函数
    public double calculateBandwidthPenalty(int[] bandwidths) {
        double bandwidthPenalty = 0.0;
        for (double bandwidth : bandwidths) {
            // 计算每个平台的带宽惩罚项
            double b_target = b_max;
            bandwidthPenalty += Math.pow((b_target - bandwidth), 2);
        }
        // 计算总的带宽限制函数，并进行归一化
        return (lambda * bandwidthPenalty) / F_bandwidth_max;  // 使用归一化的带宽惩罚

    }
}