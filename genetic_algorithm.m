clear all;clc;

% set function
f = @(var) 100*sqrt(abs(var(2) - 0.01*var(1).^2)) + 0.01*abs(var(1) + 10);

% set parameters
generations = 100;
binary_size = 6;
variables_number = 2;
population_size = 10;
crossing_bit = 6; % recomended => (binary_size*variables_number)/2;

final_size = binary_size*variables_number;
elite = 1;

genotypes = zeros(final_size, population_size);
phenotypes = zeros(variables_number, population_size);
functions = zeros(1, population_size);
probabilities = zeros(1, population_size);

% first generation
population_counter = 0;
while population_counter < population_size
    population_counter = population_counter + 1;
    genotypes_counter = 0;
    pheno = '';
    func = nan;
    phenotypes_counter = 0;
    variable_counter = 0;
    while genotypes_counter < final_size 
        genotypes_counter = genotypes_counter + 1;
        if(rand()*100 >= 50)
            gene = 1;
        else
            gene = 0;
        end
        pheno = strcat(num2str(pheno),num2str(gene));
        genotypes(genotypes_counter, population_counter) = gene;
    
        phenotypes_counter = phenotypes_counter + 1;
        if phenotypes_counter == binary_size
            phenotypes_counter = 0;
            variable_counter = variable_counter + 1;
            phenotypes(variable_counter, population_counter) = bin2dec(pheno);
            pheno = '';
        end
        
        functions(1, population_counter) = f(phenotypes(:,population_counter));
        if functions(1, elite) < functions(1, population_counter)
            elite = population_counter;
        end
    end
end

% others generations
generation = 1;
while generation < generations
    generation
    generation = generation + 1;
    
    sum_functions = sum(functions(1,:));
    population_counter = 0;
    while population_counter < population_size
        population_counter = population_counter + 1;
        probabilities(1,population_counter) = (functions(:,population_counter)/sum_functions)*100;
    end

    new_genotypes = zeros(final_size, population_size);
    new_genotypes(:, 1) = genotypes(:, elite);
    population_counter = 1;
    while population_counter < population_size
        population_counter = population_counter + 1;
        count_variables = 0;
        pair = 0;
        while true
            pair = pair + 1;
            if (probabilities(1,pair) >= (rand()*100)) && (pair ~= population_counter)
                new_genotypes(:, population_counter) = [genotypes(1:crossing_bit, pair) ; genotypes(crossing_bit+1:final_size, population_counter)];       
                break;
            end
            if pair == population_size
                pair = 0;
            end
        end
    end

    genotypes = new_genotypes;
    new_phenotypes = zeros(variables_number, population_size);
    new_functions = zeros(1, population_size);
    new_elite = 1;

    best = functions(:, elite);
    population_counter = 0;
    while population_counter < population_size
        population_counter = population_counter + 1;
        phenotypes_counter = 0;
        variable_counter = 0;
        func = nan;
        genotypes_counter = 0;
        while genotypes_counter < final_size 
            genotypes_counter = genotypes_counter + 1;
            gene = genotypes(genotypes_counter, population_counter);
            pheno = strcat(num2str(pheno),num2str(gene));

            phenotypes_counter = phenotypes_counter + 1;
            if phenotypes_counter == binary_size
                phenotypes_counter = 0;
                variable_counter = variable_counter + 1;
                new_phenotypes(variable_counter, population_counter) = bin2dec(pheno);
                pheno = '';
            end

            new_functions(1, population_counter) = f(phenotypes(:,population_counter));
            if best < new_functions(1, population_counter)
                new_elite = population_counter;
                best = new_functions(1, population_counter);
            end
        end
    end

    elite = new_elite;
    phenotypes = new_phenotypes;
    functions = new_functions;
end

disp('--------------------------------------------------')
disp('Best phenotypes: ');
disp(phenotypes(:, elite));
disp('Best genotypes: ');
disp(genotypes(:, elite));
disp('Best f(x) result: ');
disp(functions(:, elite));

x_min = -200;
x_max = 200;
y_min = -200;
y_max = 200;

x = linspace(x_min, x_max);
y = linspace(y_min, y_max);
fig = figure();
[xplot, yplot] = meshgrid(x,y);
zplot = 100*sqrt(abs(yplot - 0.01*xplot.^2)) + 0.01*abs(xplot + 10);
surf(xplot, yplot, zplot);
shading interp;
