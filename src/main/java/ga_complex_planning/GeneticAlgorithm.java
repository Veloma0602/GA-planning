package ga_complex_planning;

//import ga_complex_planning.planning_info.Info;
//import ga_complex_planning.pojo.Point;
import util.GAGraphUtil;
import util.PropertyUtil;
import java.util.*;



/**
 * @ClassName GeneticAlgorithm
 * @Description 遗传算法
 * @Author fwt
 * @Date 2024/1/25 4:35 下午
 * @Version 1.0
 **/
public class GeneticAlgorithm {

    private Properties gaComplexPro = PropertyUtil.getProperty("ga_complex");
    private Properties planningInfoPro = PropertyUtil.getProperty("planning_info");

    private Map<Double, Double> bestScoreDataMap = new HashMap<>();
    private Map<Double, Double> worstScoreDataMap = new HashMap<>();
    private Map<Double, Double> totalScoreDataMap = new HashMap<>();
    private Map<Double, Double> IterBestScoreDataMap = new HashMap<>();
    private Map<Integer, Chromosome> IterBestGeneMap = new HashMap<>();

    private double F_interference_max; // 新增：干扰的最大值
    private double F_bandwidth_max; // 新增：带宽的最大值


    public static int underZeroCount = 0;

    private static final Random r = new Random();

//
//    // 受灾地信息
//    private Map<String, Point> info = Info.getInfo();

    double wInterference = Double.valueOf(gaComplexPro.getProperty("wInterference"));
    double wBandwidth = Double.valueOf(gaComplexPro.getProperty("wBandwidth")); //权重还没设置
    double wChannel = Double.valueOf(gaComplexPro.getProperty("wChannel"));

    double ORIGIN_SCORE = Double.valueOf(gaComplexPro.getProperty("ORIGIN_SCORE"));
    // 种群
    public List<Chromosome> pop = new ArrayList<>();
    // 种群大小
    private int POP_SIZE = Integer.valueOf(gaComplexPro.getProperty("POP_SIZE"));
    // 迭代次数
    private int ITER_NUM = Integer.valueOf(gaComplexPro.getProperty("ITER_NUM"));
    // 迭代次数计数
    private int TOTAL_ITER_NUM = Integer.valueOf(gaComplexPro.getProperty("TOTAL_ITER_NUM"));
    private int iterCount = 0;
    private int TotalIterNum = 0;

    private int numOfPlatforms = Integer.valueOf(planningInfoPro.getProperty("numOfPlatforms"));

    private double f_min = Double.valueOf(planningInfoPro.getProperty("f_min"));  //频率上下限
    private double f_max = Double.valueOf(planningInfoPro.getProperty("f_max"));  //频率上下限
    private double f_disable_min = Double.valueOf(planningInfoPro.getProperty("f_disable_min"));     //禁用频率上下限
    private double f_disable_max = Double.valueOf(planningInfoPro.getProperty("f_disable_max"));
    private double b_range_min = Double.valueOf(planningInfoPro.getProperty("b_range_min"));       //带宽上下限
    private double b_range_max = Double.valueOf(planningInfoPro.getProperty("b_range_max"));
    private double b_max = Double.valueOf(planningInfoPro.getProperty("b_max"));
    private double fixedTotalPowerValue = Double.valueOf(planningInfoPro.getProperty("fixedTotalPowerValue"));    //总功率

    //记录外部循环的平均最优适应度变化
    private double[][] IterationBestValues = new double[ITER_NUM][TOTAL_ITER_NUM];
    private Map<Integer, Double> IterationAverageValuesMap = new HashMap<>();


    // 当前种群最佳适应度
    private double bestScore = Double.MIN_VALUE;
    // 当前种群最坏适应度
    private double worstScore = Double.MAX_VALUE;
    // 当前种群总适应度
    private double totalScore = 0;
    // 当前种群平均适应度
    private double averageScore = 0;
    // 当前种群最佳基因
    private Chromosome bestGene = null;
    // 当前种群最坏基因
    private Chromosome worstGene = null;

    // 变异概率
    private double MUTATION_RATE = Double.valueOf(gaComplexPro.getProperty("MUTATION_RATE"));
    // 最大变异长度
    private int MAX_MUTATION_NUM = 3;

    private CombinedPenaltyCalculator combinedPenaltyCalculator = new CombinedPenaltyCalculator();
    private double P_min = Double.MAX_VALUE;
    private double P_max = Double.MIN_VALUE;

    /**
     * 执行遗传算法函数
     * @return 执行结束后最佳的种群基因
     */
    public Chromosome conductGA() {
        double IterationBestScoreSum = 0;
        long IterTimeSum = 0;
        for(int i = 0; i < TOTAL_ITER_NUM;i++){
            long startTime = System.currentTimeMillis();
            TotalIterNum++;
            init();
            for (int j = 0; j < ITER_NUM; j++) {
                // 1、计算种群适应度
                calculatePopScore();
                print(j+1);
                // 2、交叉生成新的种群
                evolve();
                 // 3、种群变异
                 mutation();
            }
            long endTime = System.currentTimeMillis();
            System.out.println("==============================");
            System.out.println("出现负值结果个数为："+GeneticAlgorithm.underZeroCount);
            System.out.println("平均每次迭代出现负值的个数为:"+GeneticAlgorithm.underZeroCount/500);
            System.out.println("运算耗时为："+(endTime-startTime)+"ms");
            long IterTime = endTime-startTime;
            IterTimeSum += IterTime ;
            Map[] bestWorstDataSet = new Map[2];
            Map[] totalDataSet = new Map[1];
            bestWorstDataSet[0]=bestScoreDataMap;
            bestWorstDataSet[1]=worstScoreDataMap;
            totalDataSet[0]=totalScoreDataMap;
            GAGraphUtil.drawBestWorstScoreGraph(bestWorstDataSet);
            GAGraphUtil.drawTotalScoreGraph(totalDataSet);
//            return bestGene;
//            GAGraphUtil.drawCurrentPopScoreDistributeGraph(pop, iterCount);

        }

        // 初始化最大的 bestScore 为 Double.MIN_VALUE
        double maxBestScore = Double.MIN_VALUE;
        double maxBestScoreKey = Double.MIN_VALUE;
        // 遍历Map
        Iterator<Map.Entry<Double, Double>> iterator = IterBestScoreDataMap.entrySet().iterator();
        while (iterator.hasNext()) {
            Map.Entry<Double, Double> entry = iterator.next();
            double CurrentScore = entry.getValue();
            double CurrentKey = entry.getKey();
            // 如果当前的bestScore比之前的最大值大，则更新最大值
            if (CurrentScore > maxBestScore) {
                maxBestScoreKey = CurrentKey;
                maxBestScore = CurrentScore;
            }
        }

        //
        Chromosome maxBestGene = IterBestGeneMap.get((int)maxBestScoreKey);


        long IterAverTime = IterTimeSum / TOTAL_ITER_NUM;
        System.out.println("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
        System.out.println("这次运行平均运算耗时为" + IterAverTime + "ms");
        System.out.println("算法得到最佳基因的适应度为" + maxBestScore );
        System.out.println("算法得到最佳基因的频率分布为：" + Arrays.toString(maxBestGene.getFrequencies()));
        System.out.println("算法得到最佳基因的带宽分布为：" + Arrays.toString(maxBestGene.getBandwidths()));


        Map[] IterBestScoreDataSet = new Map[1];
        Map[] IterationAverageValuesSet = new Map[1];
        IterBestScoreDataSet[0] = IterBestScoreDataMap;



        //计算每次GA迭代算法平均适应度变化
        for(int n = 0; n < IterationBestValues.length; n++){
            IterationBestScoreSum = 0;
            for (int m = 0; m < IterationBestValues[n].length; m++) {
                IterationBestScoreSum +=  0.000001 * IterationBestValues[n][m];
//                System.out.println( "IterationBestScoreSum " + IterationBestScoreSum);
            }
            double average =  1000000 * (IterationBestScoreSum / IterationBestValues[n].length);
//            System.out.println("average" + average);
            IterationAverageValuesMap.put(n,average);
        }

        IterationAverageValuesSet[0] = IterationAverageValuesMap;

        GAGraphUtil.drawIterAverBestScoreGraph(IterationAverageValuesSet);

                //绘制总体GA迭代适应度变化轨迹图
        GAGraphUtil.drawIterBestScoreGraph(IterBestScoreDataSet);
        GAGraphUtil.blockUtil();

        return bestGene;

}

    /**
     * 执行遗传算法函数（收敛代数)
     *
     * @return 执行结束后最佳的种群基因

    public Chromosome conductGA() {
        double IterationBestScoreSum = 0;
        long IterTimeSum = 0;
        Chromosome maxBestGene = null;
        double maxBestScore = Double.MIN_VALUE;

        double lastBestScore = Double.MIN_VALUE; // 上一代的最佳适应度分数
        int convergenceCount = 0; // 收敛计数器
        final int convergenceThreshold = 50; // 收敛阈值，连续50代适应度变化小于1e-5视为收敛

        for (int i = 0; i < TOTAL_ITER_NUM; i++) {
            long startTime = System.currentTimeMillis();
            TotalIterNum++;
            init();
            for (int j = 0; j < ITER_NUM; j++) {
                calculatePopScore();
                print(j + 1);
                evolve();
                mutation();

                // 更新收敛判断
                if (Math.abs(bestScore - lastBestScore) < 1e-5) {
                    convergenceCount++;
                    if (convergenceCount >= convergenceThreshold) {
                        System.out.println("收敛于代数: " + (i * ITER_NUM + j + 1));
                        maxBestGene = bestGene; // 保存最佳基因
                        maxBestScore = bestScore; // 保存最佳适应度分数
                        break; // 提前终止
                    }
                } else {
                    convergenceCount = 0;
                }
                lastBestScore = bestScore;
            }
            if (maxBestGene != null) break; // 如果已经收敛，终止外部循环

            long endTime = System.currentTimeMillis();
            System.out.println("==============================");
            System.out.println("出现负值结果个数为：" + GeneticAlgorithm.underZeroCount);
            System.out.println("平均每次迭代出现负值的个数为:" + GeneticAlgorithm.underZeroCount / 500);
            System.out.println("运算耗时为：" + (endTime - startTime) + "ms");
            long IterTime = endTime - startTime;
            IterTimeSum += IterTime;
            Map[] bestWorstDataSet = new Map[2];
            Map[] totalDataSet = new Map[1];
            bestWorstDataSet[0] = bestScoreDataMap;
            bestWorstDataSet[1] = worstScoreDataMap;
            totalDataSet[0] = totalScoreDataMap;
            GAGraphUtil.drawBestWorstScoreGraph(bestWorstDataSet);
            GAGraphUtil.drawTotalScoreGraph(totalDataSet);
            // 绘图和统计数据的代码部分保留
        }

        if (maxBestGene != null) {
            System.out.println("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
            System.out.println("这次运行平均运算耗时为" + IterTimeSum / TOTAL_ITER_NUM + "ms");
            System.out.println("算法得到最佳基因的适应度为" + maxBestScore);
            System.out.println("算法得到最佳基因的频率分布为：" + Arrays.toString(maxBestGene.getFrequencies()));
            System.out.println("算法得到最佳基因的带宽分布为：" + Arrays.toString(maxBestGene.getBandwidths()));
        }

        Map[] IterBestScoreDataSet = new Map[1];
        Map[] IterationAverageValuesSet = new Map[1];
        IterBestScoreDataSet[0] = IterBestScoreDataMap;



        //计算每次GA迭代算法平均适应度变化
        for(int n = 0; n < IterationBestValues.length; n++){
            IterationBestScoreSum = 0;
            for (int m = 0; m < IterationBestValues[n].length; m++) {
                IterationBestScoreSum +=  0.000001 * IterationBestValues[n][m];
//                System.out.println( "IterationBestScoreSum " + IterationBestScoreSum);
            }
            double average =  1000000 * (IterationBestScoreSum / IterationBestValues[n].length);
//            System.out.println("average" + average);
            IterationAverageValuesMap.put(n,average);
        }

        IterationAverageValuesSet[0] = IterationAverageValuesMap;

        GAGraphUtil.drawIterAverBestScoreGraph(IterationAverageValuesSet);

                //绘制总体GA迭代适应度变化轨迹图
        GAGraphUtil.drawIterBestScoreGraph(IterBestScoreDataSet);
        GAGraphUtil.blockUtil();

        return bestGene;
    }
     */

    /**
     * 临时打印相关信息
     */
    private void print(int generation) {
        System.out.println("----------------------------------------");
        System.out.println("当前代数："+generation);
        System.out.println("当前种群最佳适应度："+bestScore);
        System.out.println("当前种群最坏适应度："+worstScore);
        System.out.println("当前种群总适应度："+totalScore);
        System.out.println("当前种群平均适应度："+ averageScore);
        System.out.println("当前种群最佳基因的频率分布："+Arrays.toString(bestGene.getFrequencies()));
        System.out.println("当前种群最佳基因的带宽分布："+Arrays.toString(bestGene.getBandwidths()));
        System.out.println("当前种群最差基因的频率分布："+Arrays.toString(worstGene.getFrequencies()));
        System.out.println("当前种群最差基因的带宽分布："+Arrays.toString(worstGene.getBandwidths()));

    }

    /**
     * 初始化种群
     */
    private void init() {
        pop.clear();
        bestScoreDataMap.clear();
        worstScoreDataMap.clear();
        totalScoreDataMap.clear();
        for (int i = 0; i < POP_SIZE ;i++) {
            Chromosome chromosome = new Chromosome(numOfPlatforms,f_min,f_max, 
            f_disable_min,f_disable_max,
            b_range_min,b_range_max,b_max);
            // 强制修正约束

            ConstraintChecker.enforceConstraints(chromosome);

            pop.add(chromosome);
        }
    }

    /**
     * 计算种群适应度
     */
    private void calculatePopScore() {
        iterCount++;
        if (pop == null || pop.size() == 0) return;
        totalScore = 0;
        averageScore = 0;
        bestScore = 0;
        worstScore = Double.MAX_VALUE;
        bestGene = null;
        worstGene = null;

        // 更新 P_min 和 P_max
        for (Chromosome chromosome : pop) {
            double rawPenalty = combinedPenaltyCalculator.calculateRawPenalty(chromosome.getFrequencies(), chromosome.getBandwidths());
            P_min = Math.min(P_min, rawPenalty);
            P_max = Math.max(P_max, rawPenalty);
        }

        for (Chromosome chromosome : pop) {
//            System.out.println(totalScore);
            calculateScore(chromosome);
            totalScore += chromosome.getScore();
            if (chromosome.getScore() > bestScore) {
                bestScore = chromosome.getScore();
                bestGene = chromosome;
            }
            if (chromosome.getScore() < worstScore) {
                worstScore = chromosome.getScore();
                worstGene = chromosome;
            }
            averageScore = totalScore / POP_SIZE;
            averageScore = averageScore > bestScore ? bestScore : averageScore;
        }
        // 将待绘制的折线图数据信息存入 map 中
        bestScoreDataMap.put((double) iterCount,bestScore);
        worstScoreDataMap.put((double) iterCount,worstScore);
        totalScoreDataMap.put((double) iterCount, totalScore);

        System.out.println("iterCount:" + iterCount);
        System.out.println("TotalIterNum:" + TotalIterNum);

        System.out.println(bestScore);

        //记录每次算法循环，种群最佳适应度在迭代过程中的平均值
        IterationBestValues[(iterCount- 1) % ITER_NUM  ][TotalIterNum - 1] = bestScore;

        //每完成一次外部循环，做一个此次算法内种群分布情况记录
        if(iterCount % ITER_NUM == 0) {
            IterBestScoreDataMap.put((double) TotalIterNum, bestScore);
            IterBestGeneMap.put(TotalIterNum,bestGene);
//            GAGraphUtil.drawCurrentPopScoreDistributeGraph(pop, iterCount);
        }

    }


    // TODO: 使用强硬的适应度计算策略 <- 1、在初始化种群的时候 2、在父代交叉的时候
    /**
     * 计算个体适应度 
     * @param chromosome
     */
     private void calculateScore(Chromosome chromosome) {
        ConstraintChecker constraintChecker = new ConstraintChecker();
        
        if (!constraintChecker.isGoodChromosome(chromosome)) {
            ConstraintChecker.enforceConstraints(chromosome);
//            throw new RuntimeException("错误：当前染色体出现错误，无法计算");
        }

        double scoreCount = 0d;


        InterferenceCalculator Icalculator = new InterferenceCalculator(chromosome);
        BandwidthConstraint Bcalculator = new BandwidthConstraint(chromosome);



         // 使用当前更新的 P_min 和 P_max 进行归一化处理
         double combinedPenalty = combinedPenaltyCalculator.calculateCombinedPenalty(chromosome.getFrequencies(), chromosome.getBandwidths(), P_min, P_max);

         double interferences = Icalculator.calculateInterferenceLevel(chromosome.getFrequencies(), chromosome.getBandwidths());
         double bandwidthPenalty = Bcalculator.calculateBandwidthPenalty(chromosome.getBandwidths());

//         System.out.println(1.0 / interferences);
//         System.out.println(bandwidthPenalty);



         // 计算适应度函数 F = w_interference * F_interference + w_bandwidth * F_bandwidth
         scoreCount = 1.0 / (wInterference * interferences + wBandwidth * (1 - bandwidthPenalty) + wChannel * combinedPenalty);
    
        // 统计出现负值分数的个数
        if (scoreCount < 0) {
            underZeroCount++;
        }
    
        // 如果适应度分数为0，要做补偿
        scoreCount += ORIGIN_SCORE;
    
        // 最终策略，防止适应度分数小于0，使轮盘赌失效
        // 如果出现大量分数为负数的情况，绘制出的图形会出现异常
        scoreCount = Math.max(0, scoreCount);
        chromosome.setScore(scoreCount);
    }
    
    

    /**
     * 轮盘赌算法  获取较优的父个体
     * @return
     */
    private Chromosome getParentChromosome() {
        // 轮盘赌选中的部分
        double slice = Math.random() * totalScore;
        //System.out.println("轮盘赌 slice:"+slice);
//        System.out.println("轮盘赌 total:"+totalScore);
        double sum = 0d;
//        int targetCount =0;
        for (Chromosome chromosome : pop) {
//            targetCount++;
            sum+=chromosome.getScore();
            // 轮盘赌选中
            if (sum>slice) {
//                System.out.println("当前命中的目标为："+targetCount);
                return chromosome;
            }
//            System.out.println("轮盘赌选中个体为:\n"+chromosome);
        }
        return pop.get(pop.size()-1);
    }

    /**
     * 交叉产生新的子代
     * 这里要注意，一定要避免迭代次数过多导致的程序阻塞
     *
     */
    private void evolve() {
        List<Chromosome> newPop = new ArrayList<>();
        while (newPop.size()<POP_SIZE) {
            Chromosome p1 = getParentChromosome();
            Chromosome p2 = getParentChromosome();
            // genetic 保证生成的子代一定合法
            List<Chromosome> children = Chromosome.genetic(p1, p2);
            ConstraintChecker constraintChecker = new ConstraintChecker();
            for (Chromosome child : children) {
                if (constraintChecker.isGoodChromosome(child)) {

                    newPop.add(child);
                }else {
//                    System.out.println("交叉时出现超出约束");
                    ConstraintChecker.enforceConstraints(child);
                    newPop.add(child);
                }
            }
        }
        // 保证新生成的子代长度与原父代长度相等
        while (newPop.size()>POP_SIZE) {
            newPop.remove(r.nextInt(newPop.size()));
        }
        pop.clear();
        pop=newPop;
    }

    /**
     * 种群变异
     * 种群变异为了保证变异后子代还是符合条件的
     * 必须对基因组序列进行两两交换
     * 
     * 保证基因组突变后还符合约束
     */
    private void mutation() {
        List<Chromosome> mutatedPop = new ArrayList<>();

        for (Chromosome chromosome : pop) {
            if (Math.random()<MUTATION_RATE) {
                Chromosome mutatedChromosome = chromosome.mutation(MAX_MUTATION_NUM,chromosome);
                ConstraintChecker constraintChecker = new ConstraintChecker();
                if (constraintChecker.isGoodChromosome(mutatedChromosome)) {

                    mutatedPop.add(mutatedChromosome);
                }else {
//                    System.out.println("变异出现越界");
                    ConstraintChecker.enforceConstraints(mutatedChromosome);  // 强制修正约束
                    mutatedPop.add(mutatedChromosome);
                }
            }else {
                mutatedPop.add(chromosome);
            }
        }

        // 保证新生成的子代长度与原父代长度相等
        while (mutatedPop.size()>POP_SIZE) {
            mutatedPop.remove(r.nextInt(mutatedPop.size()));
        }
        pop.clear();
        pop=mutatedPop;






    }
}

    // /**
    //  * 计算两点之间的路径长度
    //  * @param fromX
    //  * @param fromY
    //  * @param destX
    //  * @param destY
    //  * @return
    //  */
    // @Override
    // public double routeLength(double fromX, double fromY, double destX, double destY) {
    //     double deltaX = fromX - destX;
    //     double deltaY = fromY - destY;
    //     return Math.sqrt(deltaX*deltaX + deltaY*deltaY);
    // }

//     /**
//      * 计算基因长短
//      * 因为路径规划问题中，基因长短和车辆个数和受灾点有关系
//      * 所以需要额外开函数计算
//      * @return
//      */
//     private int calculateGeneSize() {
//         return CAR_NUM+POINT_NUM+1;
//     }

// }
