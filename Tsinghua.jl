import Pkg
using JuMP, Ipopt, CSV,DataFrames,LinearAlgebra, XLSX, IterTools, DelimitedFiles,Plots,MAT


#  Here, it is only a typical economic dispatch model

#----------------------------------Parameters----------------------------------#
P_Wmax_1=[2,1.5,1.6,1.8,1.3,0.6,2.8,3.3,3.9,4,3.3,2.9,2.7,2,0.2,3.2,5.1,3.1,1.8,2,1.3,1,2,3.8]                    # Wind_VPP1
load_1=[2.2,1.8,3,6,5.8,5.2,5.6,3.8,2.5,2.7,3,2.6,2.2,2.1,4.2,5.8,6.2,6.3,6.5,6.6,6.3,6.2,6,5.7]                  # Load_VPP1

P_Wmax_2=[4.7,5.1,4.3,4.1,3.8,3.9,4,5,5,4.8,3.9,4.3,5,5.2,5.8,5.6,1.6,0.9,5.8,4.1,3.6,3.5,3.1,3.8]                 # Wind_VPP2
load_2=[5,4,4,4.2,4.1,3.6,3.4,3.7,3.9,3.8,3.9,4,4.1,4.2,3.7,3,5.1,6.1,5.8,6.2,6.3,5.5,5,3.8]                      # Load_VPP2

P_Wmax_3=[9.3,10.1,7.2,7.5,7.9,6.4,7.1,6.9,5.6,5.4,5.2,4,3.8,3,2.8,3.2,2.5,1.1,2.1,2.9,2.7,3,4.6,5.5]             # Wind_VPP3
load_3=[4,2.1,1.1,1.1,0.7,1,1.9,3.6,3.8,4.2,5.8,5.6,5.8,5.6,5.7,6.1,8,10,9.4,8.2,6.2,5.5,4.8,2.2]                 # Load_VPP3

load_area_1=load_1+load_2
load_area_2=load_3


T=length(load_1)                      # length of market horizon

λ_BS=[0.05,0.05,0.05]                 # ES cost coefficients
P_BS_max=[0.6,0.6,1.2]                # output bounds of ES
SOC_min=[0.2,0.2,0.2]                 # min of SOC 
SOC_max=[0.9,0.9,0.9]                 # max of SOC 
SOC_initial=[0.4,0.4,0.4]             # initial SOC
BS_capacity=[1,1,2]                   # capacity of ES

b_MT=[0.9,  0.6,  0.5]                # MT cost coefficients
P_MT_initial=[2,3,4]                  # output bounds of MT
P_MT_max=[4,5,6]                      # output bounds of MT
P_MT_dn=[-2, -2.5, -3]              # climbing down bound of MT
P_MT_up=[2, 2.5, 3]                 # climbing up bound of MT

#----------------------------------Define Model----------------------------------#
model= Model()

@variable(model,0<=P_MT_1[1:T]<=P_MT_max[1])                                    #output bounds of MT
@variable(model,0<=P_MT_2[1:T]<=P_MT_max[2])
@variable(model,0<=P_MT_3[1:T]<=P_MT_max[3])

for i in 2:T
    @constraint(model, P_MT_1[i]-P_MT_1[i-1]<=P_MT_up[1])                       # bounds for MT climbing
    @constraint(model, P_MT_1[i]-P_MT_1[i-1]>=P_MT_dn[1])                       
    @constraint(model, P_MT_2[i]-P_MT_2[i-1]<=P_MT_up[2])                       
    @constraint(model, P_MT_2[i]-P_MT_2[i-1]>=P_MT_dn[2])                      
    @constraint(model, P_MT_3[i]-P_MT_3[i-1]<=P_MT_up[3])                       
    @constraint(model, P_MT_3[i]-P_MT_3[i-1]>=P_MT_dn[3])                       
end
@constraint(model, P_MT_1[1]-P_MT_initial[1]<=P_MT_up[1])
@constraint(model, P_MT_2[1]-P_MT_initial[2]<=P_MT_up[2])
@constraint(model, P_MT_3[1]-P_MT_initial[3]<=P_MT_up[3])

@variable(model, -P_BS_max[1]<=P_BS_1[1:T]<=P_BS_max[1])                    # bounds for ES output
@variable(model, -P_BS_max[2]<=P_BS_2[1:T]<=P_BS_max[2])
@variable(model, -P_BS_max[3]<=P_BS_3[1:T]<=P_BS_max[3])                    

@variable(model, SOC_min[1]<=SOC_1[1:T]<=SOC_max[1])                        # min and max of SOC 
@variable(model, SOC_min[2]<=SOC_2[1:T]<=SOC_max[2]) 
@variable(model, SOC_min[3]<=SOC_3[1:T]<=SOC_max[3]) 

@constraint(model, SOC_1[T]==SOC_initial[1])                                # the initial and end of SOC must be the same
@constraint(model, SOC_2[T]==SOC_initial[2])
@constraint(model, SOC_3[T]==SOC_initial[3])

@constraint(model, SOC_1[1]==SOC_initial[1]-P_BS_1[1]/BS_capacity[1])       # 0-1 period of SoC, initial SOC is 0.4 here
@constraint(model, SOC_2[1]==SOC_initial[2]-P_BS_2[1]/BS_capacity[2])
@constraint(model, SOC_3[1]==SOC_initial[3]-P_BS_3[1]/BS_capacity[3])

for i in 2:T
    @constraint(model, SOC_1[i]==SOC_1[i-1]-P_BS_1[i]/BS_capacity[1])           # bouns of SoC in 24 hours
    @constraint(model, SOC_2[i]==SOC_2[i-1]-P_BS_2[i]/BS_capacity[2])
    @constraint(model, SOC_3[i]==SOC_3[i-1]-P_BS_3[i]/BS_capacity[3])           
end

@variable(model,0<=P_W_1[i in 1:T]<=P_Wmax_1[i])                            # Wind turbine output upper and lower limits  
@variable(model,0<=P_W_2[i in 1:T]<=P_Wmax_2[i])
@variable(model,0<=P_W_3[i in 1:T]<=P_Wmax_3[i])                  

clearing_price=Dict()
for i in 1:T
    clearing_price[1,i]=@constraint(model, P_MT_1[i]+P_BS_1[i]+P_W_1[i] +P_MT_2[i]+P_BS_2[i]+P_W_2[i] +P_MT_3[i]+P_BS_3[i]+P_W_3[i]==load_area_1[i]+load_area_2[i])   # power balance constraint
end


C_1=b_MT[1]*sum(P_MT_1)+ λ_BS[1]*sum(P_BS_1.*P_BS_1)
C_2=b_MT[2]*sum(P_MT_2)+ λ_BS[2]*sum(P_BS_2.*P_BS_2)
C_3=b_MT[3]*sum(P_MT_3)+ λ_BS[3]*sum(P_BS_3.*P_BS_3)

@objective(model, Min, C_1+C_2+C_3)

#---------------------------solve the bi-level model-------------------
set_optimizer(model, Ipopt.Optimizer)
# set_attribute(model, "limits/gap", 0.0280)
# set_time_limit_sec(model, 700.0)
optimize!(model)



#---------------------------Output the results-------------------
P_MT_1=JuMP.value.(P_MT_1)
P_MT_2=JuMP.value.(P_MT_2)
P_MT_3=JuMP.value.(P_MT_3)
P_BS_1=JuMP.value.(P_BS_1)
P_BS_2=JuMP.value.(P_BS_2)
P_BS_3=JuMP.value.(P_BS_3)
P_W_1=JuMP.value.(P_W_1)
P_W_2=JuMP.value.(P_W_2)
P_W_3=JuMP.value.(P_W_3)
SOC_1=JuMP.value.(SOC_1)
SOC_2=JuMP.value.(SOC_2)
SOC_3=JuMP.value.(SOC_3)

lambda_1=JuMP.value.(lambda_1)
lambda_2=JuMP.value.(lambda_2)
lambda_3=JuMP.value.(lambda_3)
lambda_4=JuMP.value.(lambda_4)
lambda_5=JuMP.value.(lambda_5)
lambda_6=JuMP.value.(lambda_6)
lambda_7=JuMP.value.(lambda_7)
lambda_8=JuMP.value.(lambda_8)
lambda_9=JuMP.value.(lambda_9)
lambda_10=JuMP.value.(lambda_10)
lambda_11=JuMP.value.(lambda_11)
lambda_12=JuMP.value.(lambda_12)
lambda_13=JuMP.value.(lambda_13)
lambda_14=JuMP.value.(lambda_14)
lambda_15=JuMP.value.(lambda_15)
lambda_16=JuMP.value.(lambda_16)
lambda_17=JuMP.value.(lambda_17)
lambda_18=JuMP.value.(lambda_18)
lambda_19=JuMP.value.(lambda_19)
lambda_20=JuMP.value.(lambda_20)
lambda_21=JuMP.value.(lambda_21)
lambda_22=JuMP.value.(lambda_22)
lambda_23=JuMP.value.(lambda_23)
lambda_24=JuMP.value.(lambda_24)

lambda=dual.(values(clearing_price))
rev_comp_1=sum(lambda.*P_MT_1)
rev_comp_2=sum(lambda.*P_MT_2)
rev_comp_3=sum(lambda.*P_MT_3)

profits_1=rev_comp_1 - (b_MT[1]*sum(P_MT_1)+ λ_BS[1]*sum(P_BS_1.*P_BS_1))
profits_2=rev_comp_2 - (b_MT[2]*sum(P_MT_2)+ λ_BS[2]*sum(P_BS_2.*P_BS_2))  
profits_3=rev_comp_3 - (b_MT[3]*sum(P_MT_3)+ λ_BS[3]*sum(P_BS_3.*P_BS_3))  

energy_fee=rev_comp_2+rev_comp_3+rev_comp_1

reserve_1=zeros(T)
reserve_2=zeros(T)
reserve_3=zeros(T)
reserve_1_cons=zeros(T)
reserve_2_cons=zeros(T)
reserve_3_cons=zeros(T)
for i in 1:T
    if P_MT_1[i]>0.001
        reserve_1[i]=P_MT_max[1]-P_MT_1[i]
        if reserve_1[i]<=P_MT_up[1]
        reserve_1_cons[i]=reserve_1[i]
        else
            reserve_1_cons[i]=P_MT_up[1]
        end
    end
    if P_MT_2[i]>0.001
        reserve_2[i]=P_MT_max[2]-P_MT_2[i]  
        if reserve_2[i]<=P_MT_up[2]
            reserve_2_cons[i]=reserve_2[i]
        else
            reserve_2_cons[i]=P_MT_up[2]
        end
    end
    if P_MT_3[i]>0.001
        reserve_3[i]=P_MT_max[3]-P_MT_3[i]  
        if reserve_3[i]<=P_MT_up[3]
            reserve_3_cons[i]=reserve_3[i]
        else
            reserve_3_cons[i]=P_MT_up[3]
        end
    end

end


reserve_total=reserve_1+reserve_2+reserve_3
reserve_cons=reserve_1_cons+reserve_2_cons+reserve_3_cons

sum(reserve_total)
sum(reserve_cons)



matwrite("reserve_total_single.mat", Dict("reserve_total_single" => reserve_total))
matwrite("reserve_cons_single.mat", Dict("reserve_cons_single" => reserve_cons))

matwrite("reserve_I_comp1_single.mat", Dict("reserve_I_comp1_single" => reserve_1))
matwrite("reserve_I_comp2_single.mat", Dict("reserve_I_comp2_single" => reserve_2))
matwrite("reserve_I_comp3_single.mat", Dict("reserve_I_comp3_single" => reserve_3))

matwrite("reserve_II_comp1_single.mat", Dict("reserve_II_comp1_single" => reserve_1_cons))
matwrite("reserve_II_comp2_single.mat", Dict("reserve_II_comp2_single" => reserve_2_cons))
matwrite("reserve_II_comp3_single.mat", Dict("reserve_II_comp3_single" => reserve_3_cons))


matwrite("clearing_price_single.mat", Dict("clearing_price_single" => lambda))