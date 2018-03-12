#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>

#include <R.h>
#include <Rinternals.h>

#include <arpa/inet.h>

typedef struct {
  uint32_t poscount;
  uint8_t flag;
  float minfrac;
  uint32_t hhash;
} dump_h;

typedef struct {
  uint32_t tid;
  uint32_t p;
  uint16_t c;
  uint8_t num;
  uint16_t pos[6];
} dump_p;

size_t read_header(dump_h *h, FILE *file)
{
  size_t ret = 0;

  uint32_t poscount;
  uint8_t flag;
  uint32_t minfrac;

  ret = fread(&poscount, sizeof poscount, 1, file);
  ret = fread(&flag, sizeof flag, 1, file);
  ret = fread(&minfrac, sizeof minfrac, 1, file);

  h->poscount = ntohl(poscount);
  h->flag     = flag;
  h->minfrac  = (float)ntohl(minfrac) / 1000;

  if (h->flag) {
    uint32_t hhash;
    ret = fread(&hhash, sizeof hhash, 1, file);
    h->hhash = ntohl(hhash);
  }

  return ret;
}

size_t read_pos(dump_p *p, FILE *file, dump_h *h)
{
  size_t ret = 0;

  if (h->flag) {

    uint32_t tid;
    uint32_t pr;

    ret = fread(&tid, sizeof tid, 1, file);
    ret = fread(&pr, sizeof pr, 1, file);

    p->tid = ntohl(tid);
    p->p = ntohl(pr);
  }

  uint16_t c;
  uint8_t num;

  ret = fread(&c, sizeof c, 1, file);
  ret = fread(&num, sizeof num, 1, file);

  p->c = ntohs(c);
  p->num = num;

  if (ret == 0) // break early to avoid accessing undefined p->num
    return ret;

  uint16_t pos[6];
  memset(pos, 0, sizeof(pos));

  ret = fread(&pos, sizeof pos[0], p->num, file);

  memset(p->pos, 0, sizeof(p->pos));

  int i;
  for(i = 0; i < p->num; ++i)
    p->pos[i] = ntohs(pos[i]);

  return ret;
}


SEXP binReadR(SEXP Rin) {

  const char *in = CHAR(STRING_ELT(Rin, 0));

  FILE *file = NULL;
  file = fopen(in, "rb+");

  if (!file) {
    printf("Can't open input file: %s\n\n", in);
    return R_NilValue;
  }

  size_t ret = 0;

  dump_h dhead;
  ret = read_header( &dhead, file );


  SEXP hist = PROTECT(allocVector(VECSXP, dhead.poscount));

  int i = 0;
  while (1) {

    dump_p dpos;
    ret = read_pos( &dpos, file, &dhead);

    if (ret == 0) break;

    SET_VECTOR_ELT(hist, i, allocVector(REALSXP, dpos.num));

    int j;
    for (j = 0; j < dpos.num; ++j) {
      REAL(VECTOR_ELT(hist, i))[j] = dpos.pos[j] / (float)dpos.c;
    }

    ++i;
  }

  fclose(file);

  UNPROTECT(1);

  return hist;
}
