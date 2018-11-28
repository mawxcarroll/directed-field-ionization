private void pulseEvolution()
        {
            /*
             * Genetic algorithm happens here!
             * (1) calculate fitness score for each pulse
             * (2) erase current averaged data so next generation starts fresh
             * (3) do elitism
             * (4) do crossover
             *
             */



            /*
             * calculate fitness scores
             *  Use fraction of total signal in gate region as score
             */

            calculateFitnessScore();

            //remove last (zeroed) waveform from population before doing evolution stuff.
            population.RemoveAt(populationSize - 1);

            //update average, min, max fitness scores
            avgFitnessScore.Add(fitnessScores.Average());
            minFitnessScore.Add(fitnessScores.Min());
            maxFitnessScore.Add(fitnessScores.Max());
            avgOverlap.Add(overlapIntegrals.Average());
            minOverlap.Add(overlapIntegrals.Min());
            maxOverlap.Add(overlapIntegrals.Max());

            //reset the average data (except for the final scan)
            if (nscan < numberOfScans)
                sd.ElementAt<ScanData>(0).resetAvgData();
            //sort the pulses by fitness score
            population.Sort((x, y) => y.fitness.CompareTo(x.fitness));

            //save the previous generation first

            //if this is the first generation, create the file for writing
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


            //Console.WriteLine("----->   fitness[0] = " + population[0].fitness + ", fitness[last] = " + population[populationSize - 1].fitness);

            //print fitness scores
            /*Console.WriteLine("***********************");
            for (int i = 0; i < population.Count; i++)
            {
                Console.WriteLine(population[i].fitness);
            }
            */
            //fill up the rest of the population with new pulses
            List<FieldPulse> children = new List<FieldPulse>();
            //Add copies of the super-elites. These will be the copies that are allowed to mutate.
            for (int i = 0; i < superEliteNumber; i++)
                children.Add(population[i]);
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

                /*Console.WriteLine("tournament 1 fitness scores");
                for (int j = 0; j < tournamentSize; j++)
                    Console.WriteLine(chosen[j] + ": " + tournament[j].fitness);
                String pulse = "";
                for (int j = 0; j < parent1.pulse.Count; j++)
                {
                    pulse += Math.Round(1000*parent1.pulse[j])/1000.0 + ", ";
                    if (j % 20 == 0 && j != 0)
                        pulse += "\n";
                }
                Console.WriteLine("Parent 1:");
                Console.WriteLine(pulse);
                */
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
                if (parent1.fitness == parent2.fitness && tournamentSize > 0)
                    parent2 = tournament[1];

                /*Console.WriteLine("tournament 2 fitness scores");
                for (int j = 0; j < tournamentSize; j++)
                    Console.WriteLine(chosen[j] + ": " + tournament[j].fitness);
                pulse = "";
                for (int j = 0; j < parent2.pulse.Count; j++)
                {
                    pulse += Math.Round(1000*parent2.pulse[j])/1000.0 + ", ";
                    if (j % 20 == 0 && j != 0)
                        pulse += "\n";
                }
                Console.WriteLine("Parent 2:");
                Console.WriteLine(pulse);
                */

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
                /*
                pulse = "";
                for (int j = 0; j < child.pulse.Count; j++)
                {
                    pulse += Math.Round(1000*child.pulse[j])/1000.0 + ", ";
                    if (j % 20 == 0 && j!=0)
                        pulse += "\n";
                }
                Console.WriteLine("Child:");
                Console.WriteLine(pulse);
                */
                children.Add(child);
            }

            //add the new children to the population
            for (int i = elitismNumber; i < population.Count(); i++)
                population[i] = children[i - elitismNumber];

            //finally, mutate the population
            //VCG 1/25/17: this is where the decay of the mutation rate is set
            double progress = (double)nscan / numberOfScans;
            double currentMutationRate = mutationRate;


            /*
            //VCG 3/9/17 Changed the way the mutation rate is calculated
            double highMutationScans = 20;
            double zeroMutationScans = 5;
            if (numberOfScans <= highMutationScans + zeroMutationScans)
            {
                if (nscan <= highMutationScans)
                    currentMutationRate = mutationRate;
                else
                    currentMutationRate = 0.0;
            } // if the total number of generations is less than the combined number of high mutation generations and zero mutation generations, then just use the high mutation rate
              // for the normal number of generations (if you get that far), and then a zero mutation rate for whatever is left of the scan
            else
            {
                if (nscan <= highMutationScans)
                    currentMutationRate = mutationRate; // if we are still within the high mutation generations, the mutation rate is the user-entered value
                else if (nscan > highMutationScans && nscan <= numberOfScans - zeroMutationScans)
                    currentMutationRate = mutationRate - (nscan - highMutationScans) * (mutationRate / (numberOfScans - highMutationScans - zeroMutationScans));
                    // if we are between the high mutation generations and the zero mutation generations, the mutation rate linearly decays
                else
                    currentMutationRate = 0.0; // if we are within the zero mutation generations, mutation rate = 0
            }
            //Console.WriteLine("generation " + nscan + ", mutation rate is " + currentMutationRate);
            */

            // VCG 4/20/17: temporarily disabled the decay on the mutation rate
            /*
             if (progress < 0.1)
                 currentMutationRate = mutationRate;
             if (progress > 0.1)
             {
                 double decayRate = -0.4;//changed from -0.1 VCG 1/26/17
                 double tp = 5.0 * (progress - 0.1) / 0.9;
                 if (progress > 0.6)
                     decayRate = -0.6;//changed from -0.3 VCG 1/30/17
                 if (progress > 0.7)
                     decayRate = -0.8;//changed from -0.5 VCG 1/30/17
                 if (progress > 0.8)
                     decayRate = -1.0;
                 if (progress > 0.9)
                     decayRate = -1.5;
                 currentMutationRate = mutationRate * Math.Exp(decayRate * tp);
             }
             */

            //Console.WriteLine("mutation rate is " + currentMutationRate);
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
                population[population.Count() - 1].pulse.Add(0); // last member of population is always a zeroed waveform


        }
