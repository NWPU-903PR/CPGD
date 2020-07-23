function [scores2] = calculate_GO_scores(GO_sim_data,s1,s2);
%Function:structural scores of two modules
m=length(s1);
n=length(s2);
new_GO_matrix=GO_sim_data(s1,s2);
scores2=sum(sum(new_GO_matrix))/((m+n)*(m+n-1));



end

