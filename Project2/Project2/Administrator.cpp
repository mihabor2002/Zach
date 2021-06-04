#include "Administrator.h"
#include <string>
//211-Borovikov-Mikhail-2021
#pragma warning (disable:4996)

Administrator::Administrator() 
{}

Administrator::Administrator(const Administrator& admin) { controller = admin.controller; }

const Administrator& Administrator::operator=(const Administrator& admin) { controller = admin.controller; return *this; }

void Administrator::read_command(char *command)
{
    //Administrator(admin);
    std::ofstream log;
    log.open("log_admin.txt");
    char* c = strtok(command, " ");
    std::string com(c);
	if (com=="generate_clowd")
        	{
                log<<"zapros komandi generazii oblaka"<<std::endl;
        		if (controller.get_sost())
        		{
        			float mX,mY,sX,sY;
           			int n;
          			std::cout<<"Vvedite dispersiyu po x\n";
           			std::cin>>sX;
           			std::cout<<"Vvedite dispersiyu po y\n";
          			std::cin>>sY;
           			std::cout<<"Vvedite srednee po x\n";
           			std::cin>>mX;
         			std::cout<<"Vvedite srednee po y\n";
           			std::cin>>mY;
           			std::cout<<"Vvedite kolichestvo tochek\n";
            		std::cin>>n;
            		controller.generate_clowd(mX,mY,sX,sY,n);
            		log<<"komanda vipolnena s mX = "<<mX<<", mY = "<<mY<<", sX = "<<sX<<", sY = "<<sY<<", n = "<<n<<std::endl;
            	}
            	else
            	{
            		log<<"zapros otklonen"<<std::endl;
            		std::cout<<"komanda ne mozhet bit' vipolnena, tak kak seichas pole ne v nuzhnom sostoyanii"<<std::endl;
            	}
            }
            else if (com=="print_field")
            {
            	log<<"zapros komandi pechati polya"<<std::endl;
            	if (controller.get_sost())
            	{
            		controller.print_field();
            		log<<"komanda pechati vipolnena"<<std::endl;
            	}
            	else
            	{
          			log<<"zapros otklonen"<<std::endl;
            		std::cout<<"komanda ne mozhet bit' vipolnena, tak kak seichas pole ne v nuzhnom sostoyanii"<<std::endl;
            	}
            }
            else if (com=="voronoy")
            {
            	controller.Voronoy();
            }
            else if(com=="matrx")
            {
            	log<<"zapros komandi matrx";
            	if (controller.get_sost())
            	{
                    log<<"zapros otklonen";
                    std::cout<<"komanda ne mozhet bit' vipolnena, tak kak seichas pole ne v nuzhnom sostoyanii"<<std::endl;
            	}
            	else
            	{
            		controller.Matrx();
                    log<<"komanda vipolnena";   
            	}
            }
            else if (com=="EM")
            {
            	log<<"zapros komandi EM";
            	if (controller.get_sost())
            	{
            		log<<"zapros otklonen";
                    std::cout<<"komanda ne mozhet bit' vipolnena, tak kak seichas pole ne v nuzhnom sostoyanii"<<std::endl;
            	}
            	else
            	{
                    int k;
            		std::cin>>k;
            		controller.EM_alg(k);
            		log<<"komanda vipolnena";
            	}
            }
            else if (com=="ierarch")
            {
            	log<<"zapros komandi ierarch";
            	if (controller.get_sost())
            	{
            		log<<"zapros otklonen";
                    std::cout<<"komanda ne mozhet bit' vipolnena, tak kak seichas pole ne v nuzhnom sostoyanii"<<std::endl;            		
            	}
            	else
            	{
            		controller.Ierarch();
            		log<<"komanda vipolnena";
            	}
            }
            else if (com=="help")
            {
            	log<<"zapros komandi vivoda faila help"<<std::endl;
                std::ifstream help;
                help.open("help.txt");
                char stroka[500];
                while (!help.eof())
                {
                    help.getline(stroka,500);
                    printf("%s\n",stroka);
                }
                help.close();
                log<<"komanda vivoda faila help vipolnena"<<std::endl;
            }
            else if (com=="change_sost")
            {
            	log<<"zapros komandi izmeneniya sostoyaniya polya";
            	controller.change_sost();
            	log<<"komanda vipolnena";
            }
            else if (com=="binary")
    		{
    			log<<"zapros komandi binary";
    			if (controller.get_sost())
    			{
            		log<<"zapros otklonen"<<std::endl;
            		std::cout<<"komanda ne mozhet bit' vipolnena, tak kak seichas pole ne v nuzhnom sostoyanii"<<std::endl;
    			}
    			else
    			{
    				std::cout << "alph"<<std::endl;
    				float alph;
    				std::cin>>alph;
    				controller.Binary(alph);
    				log<<"komanda vipolnena";
    			}
    		} 
    		else if (com=="wave")
    		{
    			log<<"zapros komandi wave";
    			if (controller.get_sost())
    			{
            		log<<"zapros otklonen"<<std::endl;
            		std::cout<<"komanda ne mozhet bit' vipolnena, tak kak seichas pole ne v nuzhnom sostoyanii"<<std::endl;
                }
                else
                {
    				controller.Wave();
    				log<<"komanda vipolnena";
    			}
    		}
    		else if (com=="dbbinary")
    		{
    			log<<"zapros komandi dbbinary";
    			if (controller.get_sost())
    			{
            		log<<"zapros otklonen"<<std::endl;
            		std::cout<<"komanda ne mozhet bit' vipolnena, tak kak seichas pole ne v nuzhnom sostoyanii"<<std::endl;
                }
                else
                {
    				std::cout<<"alph"<<std::endl;
    				float alph;
    				std::cin>>alph;
                	std::cout<<"k"<<std::endl;
                	int k;
                	std::cin>>k;
    				controller.Dbbinary(alph,k);
    				log<<"komanda vipolnena";
    			}
    		}
    		else if (com=="dbwave")
    		{
    			log<<"zapros komandi dbwave";
    			if (controller.get_sost())
    			{
            		log<<"zapros otklonen"<<std::endl;
            		std::cout<<"komanda ne mozhet bit' vipolnena, tak kak seichas pole ne v nuzhnom sostoyanii"<<std::endl;
				}
				else
				{
    				controller.Dbwave();
    				log<<"komanda vipolnena";
    			}
    		}
            else if (com == "interp")
            {
            cout << "Vvedite tochku\n";
            double x,y;
            cin >> x >> y;
            Point P;
            P.x = x;
            P.y = y;
            controller.Inter(P);
 }
            else if (com == "triang")
            {
            controller.Triang();
 }
    		else if(com=="print_clasters1")
    		{
    			log<<"zapros komandi print_clasters1";
    			controller.print_clasters1();
    			log<<"komanda vipolnena";
    		}
    		else if(com=="k_sredn")
    		{
    			log<<"zapros komandi k_sredn"<<std::endl;;
    			if (controller.get_sost())
    			{
            		log<<"zapros otklonen"<<std::endl;
            		std::cout<<"komanda ne mozhet bit' vipolnena, tak kak seichas pole ne v nuzhnom sostoyanii"<<std::endl;
    			}
    			else
    			{
	    			std::cout<<"k"<<std::endl;
   		 			int k;
    				std::cin>>k;
    				controller.K_sredn(k);
    				log<<"komanda vipolnena"<<std::endl;
    			}
    		}
    		else if(com=="k_sredn_yadr")
    		{
    			log<<"zapros komandi k_sredn_yadr"<<std::endl;
    			if (controller.get_sost())
    			{
            		log<<"zapros otklonen"<<std::endl;
            		std::cout<<"komanda ne mozhet bit' vipolnena, tak kak seichas pole ne v nuzhnom sostoyanii"<<std::endl;
    			}
    			else
    			{
   		 			int k,p;
    				std::cin>>k>>p;
    				controller.K_sredn_yadr(k,p);
    				log<<"komanda vipolnena s k = "<<k<<", p = "<<p<<std::endl;
    			}
    		}
    		else if (com=="forel")
    		{
    			log<<"zapros komandi forel";
    			float R;
    			std::cin>>R;
    			controller.Forel(R);
    			log<<"komanda vipolnena";
    		}
            else if (com == "gist")
            {
                controller.gist();
            }
            else if (com=="exit")
            {
            	log<<"vihod iz programmi"<<std::endl;
            }
            else
            {
            	std::cout<<"Chto-to poshlo ne tak"<<std::endl;
            }
}