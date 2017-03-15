from cdis_pipe_utils import postgres


class MuTect(postgres.ToolTypeMixin, postgres.Base):

    __tablename__ = 'mutect2_vc_metrics'
