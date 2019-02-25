private void pulseEvolution()
        {
            /*
             * Genetic algorithm happens here!
             * (1) calculate fitness score for each pulse
             * (2) erase current averaged data so next generation starts fresh
             * (3) do elitism (we created our own custom version of elitism)
             * (4) do crossover (we used uniform crossover)
             *
             */



            /*
             * calculate fitness scores
             *  The fitness score drives the evolution of your solution. It will
             *  be idiosyncratic to your particular problem. In our case, we
             *  are trying to shape our time resolved field ionization signals
             *  to separate two overlapping signals. We tried and tested many
             *  candidate fitness scores and selected the best performer. Of
             *  course, we may not have found an optimal solution, but we really
             *  just needed one that worked.
             */

            calculateFitnessScore();

            /*
             *  remove last waveform from population before doing evolution.
             *  For the last waveform of every generation, we used an unperturbed
             *  waveform. This served as our experimental control, but we don't
             *  want to consider it in the evolution.
             */

            population.RemoveAt(populationSize - 1);

            //update average, min, max fitness scores (for display in
            // data acquisition software)
            avgFitnessScore.Add(fitnessScores.Average());
            minFitnessScore.Add(fitnessScores.Min());
            maxFitnessScore.Add(fitnessScores.Max());
            avgOverlap.Add(overlapIntegrals.Average());
            minOverlap.Add(overlapIntegrals.Min());
            maxOverlap.Add(overlapIntegrals.Max());

            //reset the average data (except for the final scan)
            //for data acquisition software display
            if (nscan < numberOfScans)
                sd.ElementAt<ScanData>(0).resetAvgData();

            /*
             *  Sort the pulses by fitness score
             *    It's not important that the fitness score be linear or
             *    in some way directly proportional to the performance. As
             *    long as the fitness score sorts the results from best to
             *    worst, it will work. In other words, only the rank of
             *    a particular solution is important.
             */

            population.Sort((x, y) => y.fitness.CompareTo(x.fitness));

            //save the previous generation first -- we write out every
            //solution for analysis

            //if this is the first generation, create the file for writing
            //this bit is particular to our data acquisition software
            if (nscan == 1)
            {
                int day = System.DateTime.Now.Day;
                int year = System.DateTime.Now.Year;
                int month = System.DateTime.Now.Month;
                popfname = saveRoot + Convert.ToString(year) + Convert.ToString(month).PadLeft(2, '0') + Convert.ToString(day).PadLeft(2, '0') + Convert.ToString(fileCounter) + ".pop";
                savePopFileStream = File.OpenWrite(popfname);
                savePopFileWriter = new BinaryWriter(savePopFileStream);
                savePopFileWriter.Write((double)population.Count());
                savePopFileWriter.Write((double)population[0].pulse.Count());
                savePopFileWriter.Write((double)geneSize);

                //for writing the best members of the final generation
                winnersfname = saveRoot + Convert.ToString(year) + Convert.ToString(month).PadLeft(2, '0') + Convert.ToString(day).PadLeft(2, '0') + Convert.ToString(fileCounter) + "_best.pop";
                saveWinnersStream = File.OpenWrite(winnersfname);
                saveWinnersWriter = new BinaryWriter(saveWinnersStream);
                //create the header
                int numToWrite = elitismNumber + superEliteNumber;
                saveWinnersWriter.Write(numToWrite); // only write the elites + the super elites
                saveWinnersWriter.Write(population[0].pulse.Count()); //write length of pulse
                saveWinnersWriter.Write(geneSize);//write size of gene
            }

            //save all members of population (except for the last)

            if (population.Count() > 0)
            {
                for (int i = 0; i < population.Count(); i++)
                    for (int j = 0; j < population[i].pulse.Count(); j++)
                        savePopFileWriter.Write(population[i].pulse[j]);
            }


            /*
             *  Generate a new population for the next generation
             */
            List<FieldPulse> children = new List<FieldPulse>();

            /*
             * Elitism and "super" Elitism
             *  Elitism is a common GA strategy. Since crossover (mating) is to
             *  some extent random and mutation is entirely random, it's possible
             *  that a new generation could be worse than the previous. To
             *  prevent this, the best members of the current population are
             *  propagated unchanged into the new population.
             *
             *  We modified this slightly. We created "super" elites that were
             *  allowed to mutate before propagating. We found this gave us a
             *  good balance between convergence and diversity. The "normal"
             *  elites were treated in the standard way and not mutated.
             */

            //Add copies of the super-elites. These will be the copies that are allowed to mutate.
            for (int i = 0; i < superEliteNumber; i++)
                children.Add(population[i]);

            /*
             * Do Crossover (mating)
             *  All pulses are eligible to mate and produce new child pulses for
             *  the next generation. The goal here is to balance diversity in
             *  your population with convergence to a good solution. We opted
             *  for "tournament selection," in which a random set of pulses
             *  are selected and then the best member of that set is chosen as
             *  one parent. This is repeated to find the second parent. The size
             *  of the tournament determines the balance between diversity and
             *  convergence (large tournaments will tend to find the highest
             *  scoring pulses).
             *
             *  We mate the pulses using "uniform crossover." For each gene, we
             *  randomly select a parent and copy that gene into the child. We
             *  have also tried "single point crossover" in which a single
             *  gene is chosen as a pivot. Genes before that pivot are copied
             *  from one parent and genes after that pivot are copied from the
             *  other parent.
             */
            for (int i = elitismNumber + superEliteNumber; i < population.Count(); i++)
            {
                //choose first parent using tournament selection
                List<FieldPulse> tournament = new List<FieldPulse>();
                List<int> chosen = new List<int>();
                int r = getRandomInt(0, population.Count() - 1);
                tournament.Add(population[r]);
                chosen.Add(r);
                for (int j = 0; j < tournamentSize; j++)
                {
                    r = getRandomInt(0, population.Count() - 1);
                    while (chosen.Contains(r))
                        r = getRandomInt(0, population.Count() - 1);
                    tournament.Add(population[getRandomInt(0, population.Count() - 1)]);
                    chosen.Add(r);
                }
                tournament.Sort((x, y) => y.fitness.CompareTo(x.fitness));
                FieldPulse parent1 = tournament[0];

                //choose second parent using tournament selection
                tournament = new List<FieldPulse>();
                chosen = new List<int>();
                r = getRandomInt(0, population.Count() - 1);
                tournament.Add(population[r]);
                chosen.Add(r);
                for (int j = 0; j < tournamentSize; j++)
                {
                    r = getRandomInt(0, population.Count() - 1);
                    while (chosen.Contains(r))
                        r = getRandomInt(0, population.Count() - 1);
                    tournament.Add(population[getRandomInt(0, population.Count() - 1)]);
                    chosen.Add(r);
                }
                tournament.Sort((x, y) => y.fitness.CompareTo(x.fitness));
                FieldPulse parent2 = tournament[0];
                //this is a bit of a hack; shouldn't treat fitness score as unique identifier...
                //(don't want identical parents!)
                if (parent1.fitness == parent2.fitness && tournamentSize > 0)
                    parent2 = tournament[1];
                //create the child using crossover: uniform or single point
                FieldPulse child = new FieldPulse();

                if (uniformCrossoverFlag)
                {
                    //do uniform crossover
                    for (int j = 0; j < parent1.pulse.Count; j++)
                    {
                        double check = getRandom(0, 1);
                        if (check < 0.5)
                            child.pulse.Add(parent1.pulse[j]);
                        else
                            child.pulse.Add(parent2.pulse[j]);
                    }
                }
                else
                {
                    //do single point crossover
                    int xpoint = getRandomInt(0, parent1.pulse.Count() - 1);
                    double check = getRandom(0, 1);
                    if (check < 0.5)
                    {
                        for (int j = 0; j < xpoint; j++)
                            child.pulse.Add(parent1.pulse[j]);
                        for (int j = xpoint; j < parent1.pulse.Count; j++)
                            child.pulse.Add(parent2.pulse[j]);
                    }
                    else
                    {
                        for (int j = 0; j < xpoint; j++)
                            child.pulse.Add(parent2.pulse[j]);
                        for (int j = xpoint; j < parent1.pulse.Count; j++)
                            child.pulse.Add(parent1.pulse[j]);
                    }

                }

                children.Add(child);
            }

            //add the new children to the population
            for (int i = elitismNumber; i < population.Count(); i++)
                population[i] = children[i - elitismNumber];

            /*
             * Finally, mutate the population
             *  We are currently using a fixed mutation rate. We have also tried
             *  a "dynamic" mutation rate that decreases as time passes.
             */
            double progress = (double)nscan / numberOfScans;
            double currentMutationRate = mutationRate;



            for (int i = superEliteNumber; i < population.Count(); i++)
            {
                for (int j = 0; j < pulseLength; j++)
                {
                    double dice = getRandom(0, 1);
                    if (dice < currentMutationRate)
                    {
                        population[i].pulse[j] = getRandom(-1 * maxSlew, maxSlew);
                    }
                }
            }

            //add back in the zeroed waveform, last member of the population

            population.Add(new FieldPulse());

            for (int j = 0; j < pulseLength; j++)
                population[population.Count() - 1].pulse.Add(0);

        }
