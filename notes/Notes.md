Exemple program C:

https://gist.github.com/ahonkela/8b23de162803d165c417

https://gist.github.com/dpryan79/b4bca03f62a02bcb6fe4


Structure d'un programme en C comptant le nombre de reads:

```C

en tete C include....


long nReads=0;

r = samopen(file);

bam1_t *b = NULL; // struct contenant 1 alignement

tant qu'on peut lire un alignement hile(sam_read1(fp, header, b) >= 0) nReads++;

samclose(r);

printf(" .....",)

return 0
``


