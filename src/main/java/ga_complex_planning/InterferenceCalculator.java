package ga_complex_planning;

import java.util.Properties;

import util.PropertyUtil;

/**
 * @ClassName InterferenceCalculator
 * @Description 干扰水平计算函数
 * @Author fwt
 * @Date 2024/1/25 4:35 下午
 * @Version 1.0
 **/

public class InterferenceCalculator {
    private double fixedTotalPowerValue; // 用于存储总功率的固定值
    private int[] frequencies; // 用于存储每个移动平台的频率

    private double F_interference_max; // 新增：干扰的最大值
    private double F_bandwidth_max; // 新增：带宽的最大值

    private double total_i = 0.0;

    private Properties gaComplexPro = PropertyUtil.getProperty("ga_complex");
    private Properties planningInfoPro = PropertyUtil.getProperty("planning_info");

    private int numOfPlatforms = Integer.valueOf(planningInfoPro.getProperty("numOfPlatforms"));
    private double b_range_min = Double.valueOf(planningInfoPro.getProperty("b_range_min"));       //带宽上下限
    private double b_range_max = Double.valueOf(planningInfoPro.getProperty("b_range_max"));


    private double[] noises = new double[numOfPlatforms];


    // 干扰水平函数中的常数
    private double alpha = Double.valueOf(gaComplexPro.getProperty("alpha"));
    //分级噪声的常数ki
    private double Ki = Double.valueOf(gaComplexPro.getProperty("Ki"));
    //电磁干扰权重
    private double E = Double.valueOf(gaComplexPro.getProperty("E"));
    //电磁干扰的频率
    private double fEMI = Double.valueOf(gaComplexPro.getProperty("fEMI"));
    //控制电磁干扰影响范围的标准差参数
    private double sigmaEMI = Double.valueOf(gaComplexPro.getProperty("sigmaEMI"));
    //大风干扰权重
    private double W = Double.valueOf(gaComplexPro.getProperty("W"));
    private double fWind = Double.valueOf(gaComplexPro.getProperty("fWind"));
    private double sigmaWind = Double.valueOf(gaComplexPro.getProperty("sigmaWind"));

    
    // 构造函数
    public InterferenceCalculator(Chromosome chromosome) {
        this.fixedTotalPowerValue = chromosome.getFixedTotalPowerValue();
        this.frequencies = chromosome.getFrequencies();
        double f_min = calculateFMin();  // 计算 f_min
        this.F_interference_max = calculateInterferenceMax(b_range_min,b_range_max, fixedTotalPowerValue, Ki, E, W, alpha, f_min);
    }

    // 计算干扰水平函数
    public double calculateInterferenceLevel(int frequencies[],int bandwidths[]) {
        for (int i = 0; i < frequencies.length; i++) {
            // 计算每个移动平台的噪声项
            double noise = calculateNoise(frequencies[i],bandwidths[i]);
            
            // Store the noise in the frequencies array
//            System.out.println(noise);

            noises[i] = noise;
            total_i += noises[i];
        }
        // 新增：归一化干扰总值
        double normalizedTotalI = alpha * total_i / F_interference_max;
        return normalizedTotalI;

    }

    //计算最大干扰值
    public double calculateInterferenceMax(double b_range_min, double b_range_max,double P_total, double Ki, double E, double W, double alpha, double f_min) {
        double NI_max = (P_total / b_range_min) + (Ki / f_min) * (1 + E + W);  // 计算每个平台的最大噪声
        double F_interference_max = alpha * NI_max * numOfPlatforms ;  // 计算总干扰的最大值
        return F_interference_max;
    }


    // 计算噪声项
    private double calculateNoise(double frequency,double bandwidth) {
        // 计算第 i 个平台的总噪声 NI_i，结合白噪声和分级噪声
        double whiteNoise = calculateWhiteNoise(bandwidth);
        double gradedNoise = calculateGradedNoise(frequency);

        // 返回总噪声
        return whiteNoise + gradedNoise;
    }

    // 计算白噪声
    private double calculateWhiteNoise(double bandwidth) {
        // 在这里实现根据频率计算白噪声的逻辑
        // 这里使用 N_{0_i}=\frac{P_{\mathrm{total_i}}}{b_i} 公式

        return  fixedTotalPowerValue / bandwidth;
    }

    // 计算分级噪声
    private double calculateGradedNoise(double frequency) {
        int platformIndex = getPlatformIndexForFrequency(frequency);
        if (platformIndex != -1) {
            // 计算 N_{f_i}，根据需要补充具体的计算逻辑
            double levelNoise = Ki / frequency * (1 + E * Math.exp(-Math.pow((frequency - fEMI), 2) / (2 * Math.pow(sigmaEMI, 2)))
                    + W * Math.exp(-Math.pow((frequency - fWind), 2) / (2 * Math.pow(sigmaWind, 2))));
            return levelNoise;
        } else {
            // 处理未找到匹配频率的情况
            return 0.0; // 或者使用其他适当的默认值
        }
    }

    // 辅助方法，根据频率获取平台索引
    private int getPlatformIndexForFrequency(double frequency) {
        for (int i = 0; i < frequencies.length; i++) {
            if (frequencies[i] == frequency) {
                return i;
            }
        }
        return -1;
    }

    //辅助归一化方法，动态初始化f_min

    private double calculateFMin() {
        // 返回 fEMI 和 fWind 中的较小值
        return Math.min(fEMI, fWind);
    }


}
