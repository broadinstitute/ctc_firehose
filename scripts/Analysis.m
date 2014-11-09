function Analysis_VA_Fig3
load ../Correlations/Correlation_and_Coverage_results.mat
load ../Correlations/SampleNames.mat
load Sample_index.mat

tumorIndex=[1:3,5:28,30:43,45:54,56:59];

for si=1:length(tumorIndex)
    currFractions=NormalizedCoverage{tumorIndex(si)}(1,:);
    currSample=SM_Names{si}(1,:);
    TotalDepth=60;
    OptimalDetection60x=OptimalDetectionSensitivity(currFractions,TotalDepth);
    TotalDepth=30;
    OptimalDetection30x=OptimalDetectionSensitivity(currFractions,TotalDepth);
    clonalFraction=1;
    Detection_vs_AvgDepth_clonal=DetectionSensitivity_vs_AvgDepth(currFractions,clonalFraction);
    clonalFraction=0.6;
    Detection_vs_AvgDepth_0_6=DetectionSensitivity_vs_AvgDepth(currFractions,clonalFraction);
    clonalFraction=0.3;
    Detection_vs_AvgDepth_0_3=DetectionSensitivity_vs_AvgDepth(currFractions,clonalFraction);
    save([currSample '_pool_analysis.mat'],'OptimalDetection30x','OptimalDetection60x','Detection_vs_AvgDepth_clonal','Detection_vs_AvgDepth_0_3','Detection_vs_AvgDepth_0_6');
end

function Detection=DetectionSensitivity_vs_noSamples(currFractions,clonalFraction)
%% Detection sensitivity for multiple libraries at different average depth per library
D=[0.25,0.5,1,2,5,10];
Fraction=currFractions(round(1./D*100));
TotalDepth=120;
individual_depths=D;
% Clonal
if clonalFraction==1
    for i=1:length(individual_depths)
        curr_depth=individual_depths(i);
        curr_detection=Fraction(i);
        clear Samples Detection_sensitivity
        MaxSamples=floor(TotalDepth/curr_depth);
        Samples=2:1:MaxSamples;
        for j=1:length(Samples)
            M=Samples(j);
            S=0;
            m=2;
            for k=0:1:m-1
                S=S+prod(M-k+1:1:M)/factorial(k)*curr_detection^(k)*(1-curr_detection)^(M-k);
            end
            Detection_sensitivity(j)=1-S;
        end
        Detection{i}=dataset(Samples',Detection_sensitivity','VarNames',{'no_samples','sensitivity'});
    %    plot(Samples*curr_depth,Detection_sensitivity,'o','MarkerEdgeColor','k','MarkerFaceColor',[1,1,1]*i/length(individual_depths));
    end
else
    cell_f=clonalFraction;
    for i=1:length(individual_depths)
        curr_depth=individual_depths(i);
        curr_detection=Fraction(i);
        clear Samples Detection_sensitivity Sensitivity
        MaxSamples=floor(TotalDepth/curr_depth);
        Samples=2:1:MaxSamples;
        for j=1:length(Samples)
                N=Samples(j);
                for M=2:1:N;
                    S=0;
                    m=2;
                    for k=0:1:m-1
                        S=S+prod(M-k+1:1:M)/factorial(k)*curr_detection^(k)*(1-curr_detection)^(M-k);
                    end
                    Sensitivity(M)=1-S;
                end
                Detection_sensitivity(j)=0;
                for M=2:1:N
                    Detection_sensitivity(j)=Detection_sensitivity(j)+exp(sum(log(N-M+1:1:N))-sum(log(1:1:M)))*cell_f^M*(1-cell_f)^(N-M)*Sensitivity(M);
                end
        end
        Detection{i}=dataset(Samples',Detection_sensitivity','VarNames',{'no_samples','sensitivity'});
        %plot(Samples*curr_depth,Detection_sensitivity,'o','MarkerEdgeColor','k','MarkerFaceColor',[1,1,1]*i/length(individual_depths));
    end
end

function Detection=DetectionSensitivity_vs_AvgDepth(currFractions,clonalFraction)
%% Detection sensitivity for subclonal variants in multiple libraries
TotalDepth=30;
D=[0.25,0.3,0.6,1,1.5,2,3,5,6,10,15];
Fraction=currFractions(round(1./D*100));
individual_depths=D;
if clonalFraction==1
% Clonal
    cell_f=1;
    for i=1:length(D)
        curr_depth=D(i);
        curr_detection=Fraction(i);
        M=floor(TotalDepth/curr_depth);
        S=0;
        m=2;
        for k=0:1:m-1
            S=S+prod(M-k+1:1:M)/factorial(k)*curr_detection^(k)*(1-curr_detection)^(M-k);
        end
        Detection_sensitivity(i)=1-S;
    end
else
    % Sub-clonal
    cell_f=clonalFraction;
    for i=1:length(individual_depths)
        curr_depth=individual_depths(i);
        curr_detection=Fraction(i);
        N=floor(TotalDepth/curr_depth);
        clear Sensitivity
        for M=2:1:N;
            S=0;
            m=2;
            for k=0:1:m-1
                S=S+prod(M-k+1:1:M)/factorial(k)*curr_detection^(k)*(1-curr_detection)^(M-k);
            end
            Sensitivity(M)=1-S;
        end
        Detection_sensitivity(i)=0;
        for M=2:1:N
            Detection_sensitivity(i)=Detection_sensitivity(i)+exp(sum(log(N-M+1:1:N))-sum(log(1:1:M)))*cell_f^M*(1-cell_f)^(N-M)*Sensitivity(M);
        end
    end
end
Detection=dataset(individual_depths',Detection_sensitivity','VarNames',{'AvgDepth','sensitivity'});

function Detection=OptimalDetectionSensitivity(currFractions,TotalDepth)
%% Optimal detection sensitivity vs. total depth
D=0.1:0.05:2;
cell_fractions=0.1:0.05:0.95;
Fraction=currFractions(round(1./D*100));
individual_depths=D;
for s=1:length(cell_fractions)
    cell_f=cell_fractions(s);
    for i=1:length(individual_depths)
        curr_depth=individual_depths(i);
        curr_detection=Fraction(i);
        N=floor(TotalDepth/curr_depth);
        clear Sensitivity
        for M=2:1:N;
            S=0;
            m=2;
            for k=0:1:m-1
                S=S+prod(M-k+1:1:M)/factorial(k)*curr_detection^(k)*(1-curr_detection)^(M-k);
            end
            Sensitivity(M)=1-S;
        end
        Detection_sensitivity(i)=0;
        for M=2:1:N
            Detection_sensitivity(i)=Detection_sensitivity(i)+exp(sum(log(N-M+1:1:N))-sum(log(1:1:M)))*cell_f^M*(1-cell_f)^(N-M)*Sensitivity(M);
        end
    end
    [val,ind]=max(Detection_sensitivity);
    optimal_depth(s)=individual_depths(ind);
    optimal_sensitivity(s)=Detection_sensitivity(ind);
end
Detection=dataset(cell_fractions',optimal_depth',optimal_sensitivity','VarNames',{'Subclonality','OptimalDepthPerSample','MaximalSensitivity'});